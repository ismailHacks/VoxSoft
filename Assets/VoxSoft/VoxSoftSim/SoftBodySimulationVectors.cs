using System.Collections;
using System.Collections.Generic;
using UnityEngine;


//Same as SoftBodySimulation but is using Vector3 instead of arrays where an index in the array is x, y, or z 
//This makes the code simpler to read buy maye a little slower according to the guy in the video, but I don't notice much difference...
public class SoftBodySimulationVectors : IGrabbable
{
	//Tetrahedralizer data structures
	private readonly TetrahedronData tetraData;
	private readonly int[] tetIds;
	private readonly int[] tetEdgeIds;

	//Same as in ball physics
	private readonly Vector3[] pos;
	private readonly Vector3[] prevPos;
	private readonly Vector3[] nodeForce;
	private readonly Vector3[] accel;
	private readonly Vector3[] accelEdge;
	private readonly Vector3[] vel;
	private readonly Vector3[] velEdge;
	private readonly Vector3[] velEdgePrev;


	//For soft body physics using tetrahedrons
	//The volume of each undeformed tetrahedron
	private readonly float[] restVolumes;
	//The length of an undeformed tetrahedron edge
	private readonly float[] restEdgeLengths;
	private readonly float[] prevEdgeLengths;
	private readonly float[] velEdgeL;


	//private readonly float[] restEdgeLengthsOriginal;

	//Inverese mass w = 1/m where m is how much mass is connected to a particle
	//If a particle is fixed we set its mass to 0
	private readonly float[] invMass;
	//Should be global so we don't have to create them a million times
	private readonly Vector3[] gradients = new Vector3[4];
	private readonly Vector3[] gradientsFacePressure = new Vector3[4];


	//The Unity mesh to display the soft body mesh
	private Mesh softBodyMesh;

	//How many vertices (particles) and tets do we have?
	private readonly int numParticles;
	private readonly int numTets;
	private readonly int numEdges;

	//Simulation settings
	private readonly Vector3 gravity = new Vector3(0f, -9.81f, 0f);
	//private readonly Vector3 gravity = new Vector3(0f, 0f, 0f);
	//3 steps is minimum or the bodies will lose their shape  
	private readonly int numSubSteps = 3;
	//To pause the simulation
	private bool simulate = true;

	//Soft body behavior settings

	//Environment collision data 
	private readonly float floorHeight = -0.01f;
	private Vector3 halfPlayGroundSize = new Vector3(5f, 8f, 5f); 

	//Grabbing with mouse to move mesh around
	//The id of the particle we grabed with mouse
	private int grabId = -1;
	//We grab a single particle and then we sit its inverted mass to 0. When we ungrab we have to reset its inverted mass to what itb was before 
	private float grabInvMass = 0f;
	//For custom raycasting
	public int[] GetMeshTriangles => tetraData.GetTetSurfaceTriIds;
	public int GetGrabId => grabId;



	public SoftBodySimulationVectors(MeshFilter meshFilter, TetrahedronData tetraData, Vector3 startPos, float meshScale = 2f)
	{
		//Tetra data structures
		this.tetraData = tetraData;

		tetIds = tetraData.GetTetIds;
		tetEdgeIds = tetraData.GetTetEdgeIds;

		numParticles = tetraData.GetNumberOfVertices;
		numTets = tetraData.GetNumberOfTetrahedrons;
		numEdges = tetraData.GetNumberOfEdges;

		//Debug.Log(numParticles + "-" + numTets+ "-" + numEdges);

		//Init the arrays 
		//Has to be done in the constructor because readonly
		pos = new Vector3[numParticles];
		prevPos = new Vector3[numParticles];
		nodeForce = new Vector3[numParticles];
		accel = new Vector3[numParticles];
		accelEdge = new Vector3[numParticles];
		vel = new Vector3[numParticles];
		velEdge = new Vector3[numParticles];
		velEdgePrev = new Vector3[numParticles];
		invMass = new float[numParticles];

		restVolumes = new float[numTets];
		restEdgeLengths = new float[numEdges];
		prevEdgeLengths = new float[numEdges];
		velEdgeL = new float[numEdges];

		//restEdgeLengthsOriginal = new float[numEdges];

		//Fill the arrays
		FillArrays(meshScale);

		//Move the mesh to its start position
		Translate(startPos);

		//Init the mesh
		InitMesh(meshFilter, tetraData);
	}

	//Fill the data structures needed or soft body physics
	private void FillArrays(float meshScale)
	{
		//[x0, y0, z0, x1, y1, z1, ...]
		float[] flatVerts = tetraData.GetVerts;

		//Particle position
		for (int i = 0; i < flatVerts.Length; i += 3)
		{
			float x = flatVerts[i + 0];
			float y = flatVerts[i + 1];
			float z = flatVerts[i + 2];

			pos[i / 3] = new Vector3(x, y, z) * meshScale;
		}

		//Particle previous position
		//Not needed because is already set to 0s

		//Particle velocity
		//Not needed because is already set to 0s

		//Rest volume
		for (int i = 0; i < numTets; i++)
		{
			restVolumes[i] = GetTetVolume(i);
		}

		//Inverse mass (1/w)
		for (int i = 0; i < numTets; i++)
		{
			float vol = restVolumes[i];

			//The mass connected to a particle in a tetra is roughly volume / 4
			float pInvMass = vol > 0f ? 1f / (vol / 4f) : 0f;

			invMass[tetIds[4 * i + 0]] += pInvMass;
			invMass[tetIds[4 * i + 1]] += pInvMass;
			invMass[tetIds[4 * i + 2]] += pInvMass;
			invMass[tetIds[4 * i + 3]] += pInvMass;
		}

		//Rest edge length
		for (int i = 0; i < restEdgeLengths.Length; i++)
		{
			int id0 = tetEdgeIds[2 * i + 0];
			int id1 = tetEdgeIds[2 * i + 1];

			restEdgeLengths[i] = Vector3.Magnitude(pos[id0] - pos[id1]);
			prevEdgeLengths[i] = Vector3.Magnitude(pos[id0] - pos[id1]);

			//restEdgeLengthsOriginal[i] = Vector3.Magnitude(pos[id0] - pos[id1]);
		}
	}

	public void MyFixedUpdate(float springConstant, float damperConstant, float volConstant)
	{
		if (!simulate)
		{
			return;
        }

		float dt = Time.fixedDeltaTime;

		//ShrinkWalls(dt);

		Simulate(dt, springConstant, damperConstant, volConstant);
	}

	public void MyUpdate()
	{
		simulate = true;
	
		//Launch the mesh upwards when pressing space
		if (Input.GetKey(KeyCode.Space))
		{
			Yeet();
		}

		//Make the mesh flat when holding right mouse 
		if (Input.GetMouseButton(1))
		{
			Squeeze();

			simulate = false;
		}

		if (simulate)
		{
			//Update the visual mesh
			UpdateMesh();
		}
		
		for (int i = 0; i < numTets; i++)
		{	
			for (int j = 0; j < 4; j++)
			{
                //The 3 opposite vertices ids
                int id0 = tetIds[4 * i + TetrahedronData.volIdOrder[j][0]];
                int id1 = tetIds[4 * i + TetrahedronData.volIdOrder[j][1]];
                int id2 = tetIds[4 * i + TetrahedronData.volIdOrder[j][2]];
			
				Debug.DrawRay(pos[id0], gradientsFacePressure[j].normalized, Color.red);
				Debug.DrawRay(pos[id1], gradientsFacePressure[j].normalized, Color.green);
				Debug.DrawRay(pos[id2], gradientsFacePressure[j].normalized, Color.blue);
			}
		}
	}

	public Mesh MyOnDestroy()
	{
		return softBodyMesh;
	}


	void Simulate(float dt, float springConstant, float damperConstant, float volConstant)
	{
		float sdt = dt / numSubSteps;

		for (int step = 0; step < numSubSteps; step++)
		{		
			PreSolve(sdt, gravity);
			//Debug.Log(" 1. " +pos[0]+ " | " + pos[1] + " | "+ pos[2] + " | "+ pos[3] + " | ");
			SolveConstraints(sdt, springConstant, damperConstant, volConstant);

			PosMove(dt, 0.01f);

			HandleEnvironmentCollision();

			PostSolve(dt);
		}
	}
 
	void PreSolve(float dt, Vector3 gravity)
	{
		//For each particle
		for (int i = 0; i < numParticles; i++)
		{
			//This means the particle is fixed, so don't simulate it
			if (invMass[i] == 0f)
			{
				continue;
			}

			//Save old pos
			prevPos[i] = pos[i];
			nodeForce[i] = gravity*(1/invMass[i]);
			//Update vel
			//vel[i] += dt * gravity;
			//Update pos
			//pos[i] += dt * vel[i];
		}
	}

	private void SolveConstraints(float dt, float springConstant, float damperConstant, float volConstant)
	{
		SolveEdges(springConstant, damperConstant, dt);
		//SolveVolumes(volCompliance, dt);
		//SolvePressure(dt, volScale);
	}

	private void SolveEdges(float springConstant, float damperConstant, float dt)
	{
		for (int i = 0; i < numEdges; i++)
		{
			//2 vertices per edge in the data structure, so multiply by 2 to get the correct vertex index
			int id0 = tetEdgeIds[2 * i + 0];
			int id1 = tetEdgeIds[2 * i + 1];

			float w0 = invMass[id0];
			float w1 = invMass[id1];

			float wTot = w0 + w1;
			
			//This edge is fixed so dont simulate
			if (wTot == 0f)
			{
				continue;
			}

			//The current length of the edge l

			//x0-x1
			//The result is stored in grads array
			Vector3 id0_minus_id1 = pos[id0] - pos[id1];

			//sqrMargnitude(x0-x1)
			float l = Vector3.Magnitude(id0_minus_id1);

			//If they are at the same pos we get a division by 0 later so ignore
			if (l == 0f)
			{
				continue;
			}

			Vector3 gradC = id0_minus_id1 / l;
			
			float x = l - restEdgeLengths[i];
			velEdgeL[id0] += (prevEdgeLengths[i]-x)*dt;
			//Debug.Log(damperConstant + " - " + velEdgeL[0] + " - " + (velEdgeL[id0]*damperConstant*gradC));
			prevEdgeLengths[id0] = x;

			//Velocity change not working due to not just edge damping taking place
			//nodeForce[id0] = nodeForce[id0] - (springConstant*C*gradC) - (vel[id0]*damperConstant);
			//nodeForce[id1] = nodeForce[id1] + (springConstant*C*gradC) - (vel[id0]*damperConstant);
			/*accelEdge[id0] = (nodeForce[id0] - (springConstant*x*gradC))*invMass[id0];
			velEdge[id0] += accelEdge[id0]*dt;

			accelEdge[id1] = (nodeForce[id1] + (springConstant*x*gradC))*invMass[id1];
			velEdge[id1] += accelEdge[id1]*dt;*/

			nodeForce[id0] = nodeForce[id0] - (springConstant*x*gradC) - (velEdgeL[id0]*damperConstant*gradC);
			nodeForce[id1] = nodeForce[id1] + (springConstant*x*gradC) + (velEdgeL[id0]*damperConstant*gradC);
		}
	}

	private void SolveVolumes(float volConstant, float dt)
	{
		//For each tetra
		for (int i = 0; i < numTets; i++)
		{
			float wTimesGrad = 0f;
		
			//Foreach vertex in the tetra
			for (int j = 0; j < 4; j++)
			{
				int idThis = tetIds[4 * i + j];

                //The 3 opposite vertices ids
                int id0 = tetIds[4 * i + TetrahedronData.volIdOrder[j][0]];
                int id1 = tetIds[4 * i + TetrahedronData.volIdOrder[j][1]];
                int id2 = tetIds[4 * i + TetrahedronData.volIdOrder[j][2]];

                //(x4 - x2)
                Vector3 id1_minus_id0 = pos[id1] - pos[id0];
				//(x3 - x2)
				Vector3 id2_minus_id0 = pos[id2] - pos[id0];
				//(x4 - x2)x(x3 - x2)
				Vector3 cross = Vector3.Cross(id1_minus_id0, id2_minus_id0);

				//Dividing by 6 due to the volume equation for a Tetrahedron
				//There may be something wrong here as the gradient is not normalized or scaled to account for t
				//Vector3 gradC = cross * (1f / 6f);
				Vector3 gradC = cross;

				gradients[j] = gradC;

				//w1 * |grad_C1|^2
				wTimesGrad += invMass[idThis] * Vector3.SqrMagnitude(gradC);
			}

			//All vertices are fixed so dont simulate
			if (wTimesGrad == 0f) 
			{;
				continue;
			}

			float vol = GetTetVolume(i);
			float restVol;
			restVol = restVolumes[i];

			float C = 6*(vol - restVol);
			//float lambda = -C / (wTimesGrad + alpha);
			
            //Move each vertex
            for (int j = 0; j < 4; j++)
            {
                int id = tetIds[4 * i + j];

				//Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
				//pos[id] += (lambda * invMass[id] * gradients[j]) + (volScale * gradients[j]); //Added (volScale * gradients[j]) to be able to control volume increase
				//pos[id] += lambda * invMass[id] * gradients[j]; //Added (volScale * gradients[j]) to be able to control volume increase
			}
		}
	}

	private void SolvePressure(float dt, float presScale)
	{
		//1. Find particles in a single triangle face
		//2. Find the normal to that triangle
		//3. Find the force acting on the face and the mass of each particle
		//4. Calculate acceleration and add to that face
		/*float oneOverdt = 1f / dt;
		for (int i = 0; i < numParticles; i++)
		{
			prevPos[i] = pos[i];
		}*/
		
		//Debug.Log("1. = " + vel[1].y*dt);
		
		//1 and 2
		for (int i = 0; i < numTets; i++)
		{
			Vector3[] posDeltas = new Vector3[4];
			for (int k = 0; k < 4; k++)
			{
				posDeltas[k] = new Vector3(0,0,0);
			}

			//Foreach vertex in a plane
			for (int j = 0; j < 4; j++)
			{
				int idThis = tetIds[4 * i + j];

                //The 3 opposite vertices ids
                int id0 = tetIds[4 * i + TetrahedronData.volIdOrder[j][0]];
                int id1 = tetIds[4 * i + TetrahedronData.volIdOrder[j][1]];
                int id2 = tetIds[4 * i + TetrahedronData.volIdOrder[j][2]];

                Vector3 id1_minus_id0 = pos[id1] - pos[id0];
				Vector3 id2_minus_id0 = pos[id2] - pos[id0];
				Vector3 cross = Vector3.Cross(id1_minus_id0, id2_minus_id0);

				float faceArea = cross.magnitude*0.5f;

				gradientsFacePressure[j] = cross;

				/*int[] ids = { id0, id1, id2 };

				foreach (int id in ids)
				{
					if (invMass[id] == 0)
					{
						velF[id] = new Vector3(0, 0, 0);
					}
					else
					{
						Vector3 accel = (presScale * faceArea / 3 * gradientsFacePressure[j].normalized)/invMass[id];
						//velF[id] += accel * dt;
						//Debug.Log(vel[id].x + " - " + vel[id].y + " - " + vel[id].z);
						vel[id] = vel[id]*dt + 0.5f*accel*dt*dt;
						//velF[id] = (presScale * faceArea / 3 * gradientsFacePressure[j].normalized) * dt;
						//Debug.Log(gradientsFacePressure[j]);
					}
				}*/
				
				Vector3 accelid0 = (presScale * faceArea / 3 * gradientsFacePressure[j])/invMass[id0];
				Vector3 accelid1 = (presScale * faceArea / 3 * gradientsFacePressure[j])/invMass[id1];
				Vector3 accelid2 = (presScale * faceArea / 3 * gradientsFacePressure[j])/invMass[id2];

				vel[id0] = vel[id0]*dt + accelid0*dt*dt;
				vel[id1] = vel[id1]*dt + accelid1*dt*dt;
				vel[id2] = vel[id2]*dt + accelid2*dt*dt;

				posDeltas[id0] += vel[id0];
				posDeltas[id1] += vel[id1];
				posDeltas[id2] += vel[id2];
			}

			
			for (int j = 0; j < 4; j++)
            {
                int id = tetIds[4 * i + j];
				pos[id] += posDeltas[id];
				posDeltas[id] = new Vector3(0,0,0);
			}
		}

		//Debug.Log("2. = " + vel[1].y);
		
		/*for (int i = 0; i < numParticles; i++)
		{
			if (invMass[i] == 0f)
			{
				continue;
			}
			//v = (x - xPrev) / dt
			velF[i] = (pos[i] - prevPos[i]) * oneOverdt;
		}*/
	}

	private void PosMove(float dt, float damperConstant)
	{
		float oneOverdt = 1f / dt;
		for (int i = 0; i < numParticles; i++)
		{
			//nodeForce[i] -= vel[i]*damperConstant;
			accel[i] = nodeForce[i]*invMass[i];
			vel[i] += accel[i]*dt;
			pos[i] += vel[i]*dt;
		}
	}


	//Environment collision handling
	private void HandleEnvironmentCollision()
	{
		for (int i = 0; i < numParticles; i++)
		{
			EnvironmentCollision(i);
		}
	}

	//Fix velocity
	private void PostSolve(float dt)
	{
		float oneOverdt = 1f / dt;
		//For each particle
		for (int i = 0; i < numParticles; i++)
		{
			if (invMass[i] == 0f)
			{
				continue;
			}

			//v = (x - xPrev) / dt
			vel[i] = (pos[i] - prevPos[i]) * oneOverdt;
		}
	}

	//Collision with invisible walls and floor
	private void EnvironmentCollision(int i)
	{
		//Floor collision
		float x = pos[i].x;
		float y = pos[i].y;
		float z = pos[i].z;

		//X
		if (x < -halfPlayGroundSize.x)
		{
			pos[i] = prevPos[i];
			pos[i].x = -halfPlayGroundSize.x;
		}
		else if (x > halfPlayGroundSize.x)
		{
			pos[i] = prevPos[i];
			pos[i].x = halfPlayGroundSize.x;
		}

		//Y
		if (y < floorHeight)
		{
			//Set the pos to previous pos
			pos[i] = prevPos[i];
			//But the y of the previous pos should be at the ground
			pos[i].y = floorHeight;
		}
		else if (y > halfPlayGroundSize.y)
		{
			pos[i] = prevPos[i];
			pos[i].y = halfPlayGroundSize.y;
		}

		//Z
		if (z < -halfPlayGroundSize.z)
		{
			pos[i] = prevPos[i];
			pos[i].z = -halfPlayGroundSize.z;
		}
		else if (z > halfPlayGroundSize.z)
		{
			pos[i] = prevPos[i];
			pos[i].z = halfPlayGroundSize.z;
		}
	}

	//
	// Unity mesh 
	//

	//Init the mesh when the simulation is started
	private void InitMesh(MeshFilter meshFilter, TetrahedronData tetraData)
	{
		Mesh mesh = new();
		//Debug.Log(pos[0]);
		mesh.vertices = pos;
		mesh.triangles = tetraData.GetTetSurfaceTriIds;

		mesh.RecalculateBounds();
		mesh.RecalculateNormals();

		meshFilter.sharedMesh = mesh;

		softBodyMesh = meshFilter.sharedMesh;
		softBodyMesh.MarkDynamic();
	}

	//Update the mesh with new vertex positions
	private void UpdateMesh()
	{
		softBodyMesh.vertices = pos;

		softBodyMesh.RecalculateBounds();
		softBodyMesh.RecalculateNormals();
	}

	//
	// Help methods
	//

	//Move all vertices a distance of (x, y, z)
	private void Translate(Vector3 moveDist)
	{
		for (int i = 0; i < numParticles; i++)
		{
			pos[i] += moveDist;
			prevPos[i] += moveDist;
		}
	}

	//Calculate the volume of a tetrahedron
	private float GetTetVolume(int nr)
	{
		//The 4 vertices belonging to this tetra 
		int id0 = tetIds[4 * nr + 0];
		int id1 = tetIds[4 * nr + 1];
		int id2 = tetIds[4 * nr + 2];
		int id3 = tetIds[4 * nr + 3];

		Vector3 a = pos[id0];
		Vector3 b = pos[id1];
		Vector3 c = pos[id2];
		Vector3 d = pos[id3];

		float volume = Tetrahedron.Volume(a, b, c, d);

		return volume;
	}

	//
	// Mesh user interactions
	//

	//Shrink walls
	private void ShrinkWalls(float dt)
	{
		//Shrink walls
		float wallSpeed = 0.4f;

		halfPlayGroundSize.x -= wallSpeed * dt;
		halfPlayGroundSize.z -= wallSpeed * dt;

		float minWallSize = 0.2f;

		halfPlayGroundSize.x = Mathf.Clamp(halfPlayGroundSize.x, minWallSize, 100f);
		halfPlayGroundSize.z = Mathf.Clamp(halfPlayGroundSize.z, minWallSize, 100f);
	}

	//Yeet the mesh upwards
	private void Yeet()
	{
		Translate(new Vector3(0f, 0.2f, 0f));
	}

	//Squash the mesh so it becomes flat against the ground
	void Squeeze()
	{
		for (int i = 0; i < numParticles; i++)
		{
			//Set y coordinate to slightly above floor height
			pos[i].y = floorHeight + 0.01f;
		}

		UpdateMesh();
	}

	//Input pos is the pos in a triangle we get when doing ray-triangle intersection
	public void StartGrab(Vector3 triangleIntersectionPos)
	{
		//Find the closest vertex to the pos on a triangle in the mesh
		float minD2 = float.MaxValue;
		
		grabId = -1;
		
		for (int i = 0; i < numParticles; i++)
		{
			float d2 = Vector3.SqrMagnitude(triangleIntersectionPos - pos[i]);
			
			if (d2 < minD2)
			{
				minD2 = d2;
				grabId = i;
			}
		}

		//We have found a vertex
		if (grabId >= 0)
		{
			//Save the current innverted mass
			grabInvMass = invMass[grabId];
			
			//Set the inverted mass to 0 to mark it as fixed
			invMass[grabId] = 0f;

			//Set the position of the vertex to the position where the ray hit the triangle
			pos[grabId] = triangleIntersectionPos;
		}
	}

	public void MoveGrabbed(Vector3 newPos)
	{
		if (grabId >= 0)
		{
			pos[grabId] = newPos;
		}
	}

	public void EndGrab(Vector3 newPos, Vector3 newParticleVel)
	{
		if (grabId >= 0)
		{
			//Set the mass to whatever mass it was before we grabbed it
			invMass[grabId] = grabInvMass;

			vel[grabId] = newParticleVel;
		}

		grabId = -1;
	}

	public void IsRayHittingBody(Ray ray, out CustomHit hit)
	{
		//Mesh data
		Vector3[] vertices = pos;

		int[] triangles = GetMeshTriangles;

		//Find if the ray hit a triangle in the mesh
		Intersections.IsRayHittingMesh(ray, vertices, triangles, out hit);
	}

	public Vector3 GetGrabbedPos()
	{
		return pos[grabId];
    }
}

