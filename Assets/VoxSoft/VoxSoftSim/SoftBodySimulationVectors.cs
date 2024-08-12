using System.Collections;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;


//Same as SoftBodySimulation but is using Vector3 instead of arrays where an index in the array is x, y, or z 
public class SoftBodySimulationVectors : IGrabbable
{
	//Tetrahedralizer data structures
	private readonly TetrahedronData tetraData;
	public voxelTet vox;
	private readonly int[] tetIds;
	private readonly int[] tetEdgeIds;
	private readonly float density = 1000; //kg/m^3
	private float sensitivity = 0.00001f;
	private float startingVerticalDisplacement = 0.15f;
	private float fitnessExponential = -10f;
	public bool converged = false;

	public static int[] cubeFloorVoxels = new int[] {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
	public static int[] cubeCompressionVoxels = new int[25];



	private readonly Vector3[] pos;
	private readonly Vector3[] prevPos;
	private readonly Vector3[] vel;

	private readonly Vector3[] stabPos;
	private readonly Vector3[] stabPrevPos;
	private readonly Vector3[] stabVel;

	//For soft body physics using tetrahedrons
	//The volume of each undeformed tetrahedron
	private readonly float[] restVolumes;
	//The length of an undeformed tetrahedron edge
	private readonly float[] restEdgeLengths;
	private readonly float[] restEdgeLengthsOriginal;

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
	//To pause the simulation
	private bool simulate = true;
	//Environment collision data 
	private readonly float floorHeight = 0f;
	private Vector3 halfPlayGroundSize = new Vector3(5f, 8f, 5f); 

	//Grabbing with mouse to move mesh around
	//The id of the particle we grabed with mouse
	private int grabId = -1;
	//We grab a single particle and then we sit its inverted mass to 0. When we ungrab we have to reset its inverted mass to what itb was before 
	private float grabInvMass = 0f;
	//For custom raycasting
	public int[] GetMeshTriangles => tetraData.GetTetSurfaceTriIds;
	public int GetGrabId => grabId;



	public SoftBodySimulationVectors(MeshFilter meshFilter, TetrahedronData tetraData, Vector3 startPos)
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
		vel = new Vector3[numParticles];

		stabPos = new Vector3[numParticles];
		stabPrevPos = new Vector3[numParticles];
		stabVel = new Vector3[numParticles];

		invMass = new float[numParticles];

		restVolumes = new float[numTets];
		restEdgeLengths = new float[numEdges];
		restEdgeLengthsOriginal = new float[numEdges];

		//Fill the arrays
		FillArrays();

		//Move the mesh to its start position
		Translate(startPos);

		//Init the mesh
		InitMesh(meshFilter, tetraData);
	}

	//Fill the data structures needed or soft body physics
	private void FillArrays()
	{
		//[x0, y0, z0, x1, y1, z1, ...]
		float[] flatVerts = tetraData.GetVerts;

		//Particle position
		for (int i = 0; i < flatVerts.Length; i += 3)
		{
			float x = flatVerts[i + 0];
			float y = flatVerts[i + 1];
			float z = flatVerts[i + 2];

			pos[i / 3] = new Vector3(x, y, z);
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

		float totalMass = 0f;

		//Inverse mass (1/w)
		for (int i = 0; i < numTets; i++)
		{
			float vol = restVolumes[i];
			float tetMass = vol*density;
			totalMass += tetMass;

			//The mass connected to a particle in a tetra is roughly volume / 4
			//float pInvMass = vol > 0f ? 1f / (tetMass/4) : 0f;

			/*invMass[tetIds[4 * i + 0]] += pInvMass;
			invMass[tetIds[4 * i + 1]] += pInvMass;
			invMass[tetIds[4 * i + 2]] += pInvMass;
			invMass[tetIds[4 * i + 3]] += pInvMass;*/
		}
		totalMass = totalMass/2;
		float pMass = totalMass/numParticles;
		Debug.Log("Num Vertices = " + numParticles);
		Debug.Log("Mass = " + totalMass);

		for (int i = 0; i < numParticles; i++)
		{
			invMass[i] = 1/pMass;
		}

		for (int i = 0; i < cubeFloorVoxels.Length; i++)
		{
			cubeCompressionVoxels[i] = cubeFloorVoxels[i]+100;
		}
		
		
		//Rest edge length
		for (int i = 0; i < restEdgeLengths.Length; i++)
		{
			int id0 = tetEdgeIds[2 * i + 0];
			int id1 = tetEdgeIds[2 * i + 1];

			restEdgeLengths[i] = Vector3.Magnitude(pos[id0] - pos[id1]);
			restEdgeLengthsOriginal[i] = Vector3.Magnitude(pos[id0] - pos[id1]);
		}
	}

	public void MyFixedUpdate(int numSubSteps, float edgeCompliance, float volCompliance, float dampingCoefficient, float pressure)
	{
		if (!simulate)
		{
			return;
        }

		float dt = Time.fixedDeltaTime;

		//ShrinkWalls(dt);
		//convergenceSetup();
		Simulate(dt, numSubSteps, edgeCompliance, volCompliance, dampingCoefficient, pressure);
		//converged = convergenceDetect(sensitivity, dt);
		//fitnessCalculate();
		//Debug.Log(converged1);
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
	}

	public Mesh MyOnDestroy()
	{
		return softBodyMesh;
	}

	//
	// Simulation
	//

	//Main soft body simulation loop
	void Simulate(float dt, int numSubSteps, float edgeCompliance, float volCompliance, float dampingCoefficient,  float pressure)
	{
		float sdt = dt / numSubSteps;

		for (int step = 0; step < numSubSteps; step++)
		{	
			PreSolve(sdt, gravity);
			SolveConstraints(sdt, edgeCompliance, volCompliance, dampingCoefficient, pressure);
			HandleEnvironmentCollision();
			PostSolve(sdt);
		}
		debugLog();
	}

	//Move the particles and handle environment collision
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
			prevPos[i] = pos[i];
			vel[i] += dt * gravity;
		}
	}

	//Handle the soft body physics
	private void SolveConstraints(float dt, float edgeCompliance, float volCompliance, float dampingCoefficient, float pressure)
	{
		//Constraints
		//Enforce constraints by moving each vertex: x = x + deltaX
		//- Correction vector: deltaX = lambda * w * gradC
		//- Inverse mass: w
		//- lambda = -C / (w1 * |grad_C1|^2 + w2 * |grad_C2|^2 + ... + wn * |grad_C|^2 + (alpha / dt^2)) where 1, 2, ... n is the number of participating particles in the constraint.
		//		- n = 2 if we have an edge, n = 4 if we have a tetra
		//		- |grad_C1|^2 is the squared length
		//		- (alpha / dt^2) is what makes the costraint soft. Remove it and you get a hard constraint
		//- Compliance (inverse stiffness): alpha

		//SolvePressureForce(dt, pressure, voxEnd, voxelTet.voxelTop);
		//SolveExternalVoxelPressureForce(dt, pressure);
		lockFaces(cubeFloorVoxels, voxelTet.voxelPositiveY);
		SolvePressureForce(dt, pressure, cubeCompressionVoxels, voxelTet.voxelNegativeY);

		forceMove(dt, dampingCoefficient);
		SolveEdges(dt, edgeCompliance);
		SolveVolumes(dt, volCompliance);
	}

	//Solve distance constraint
	private void SolveEdges(float dt, float edgeCompliance)
	{
		float alpha = edgeCompliance / (dt * dt);

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

			//(xo-x1) * (1/|x0-x1|) = gradC
			Vector3 gradC = id0_minus_id1 / l;
			float l_rest;
			
			l_rest = restEdgeLengths[i];
			//l_rest = restEdgeLengths[i]*Mathf.Pow(volScale, 1/3f);
			float C = l - l_rest;

			//lambda because |grad_Cn|^2 = 1 because if we move a particle 1 unit, the distance between the particles also grows with 1 unit, and w = w0 + w1
			float lambda = -C / (wTot + alpha);
			//float lambda = -C / (wTot);

			
			//Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
			pos[id0] += lambda * w0 * gradC;
			pos[id1] += -lambda * w1 * gradC;
		}
	}

	//Solve volume constraint
	private void SolveVolumes(float dt, float volumeCompliance)
	{
		float alpha = volumeCompliance / (dt * dt);

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
			float lambda = -C / (wTimesGrad + alpha);
			
            //Move each vertex
            for (int j = 0; j < 4; j++)
            {
                int id = tetIds[4 * i + j];

				//Move the vertices x = x + deltaX where deltaX = lambda * w * gradC
				//pos[id] += (lambda * invMass[id] * gradients[j]) + (volScale * gradients[j]); //Added (volScale * gradients[j]) to be able to control volume increase
				pos[id] += lambda * invMass[id] * gradients[j]; //Added (volScale * gradients[j]) to be able to control volume increase
			}
		}
	}

	//Used to lock specific faces in space.
	private void lockFaces(int[] voxIDs, int[] face)
	{
        for (int i = 0; i < voxIDs.Length; i++)
        {
			for (int j = 0; j < 4; j++)
            {
				int[] vertexMapping = tetraData.GetVertexMapping;
                // The id's of all particles on the face
				// Need to fix for voxID's and TetID's
                invMass[vertexMapping[8 * voxIDs[i] + face[0]]] = 0f;
                invMass[vertexMapping[8 * voxIDs[i] + face[1]]] = 0f;
                invMass[vertexMapping[8 * voxIDs[i] + face[2]]] = 0f;
                invMass[vertexMapping[8 * voxIDs[i] + face[3]]] = 0f;
			}
		}
	}

	//Damping coefficient currently implemented incorrectly - as it damps gravity as well. Just here for simulation stability
	//When dynamic effects are being looked at this will need to be correct.
    private void forceMove(float dt, float dampingCoefficient)
    {
        // Update positions based on velocity
        for (int i = 0; i < numParticles; i++)
        {
            if (invMass[i] != 0)
            {
                pos[i] += vel[i]*dampingCoefficient * dt;
            }
        }
		//Debug.DrawRay(pos[10], gravity, Color.blue);
    }

	//Update the velocity after the constrain has been handled
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
    
	//Environment collision handling
    private void HandleEnvironmentCollision()
	{
		for (int i = 0; i < numParticles; i++)
		{
			EnvironmentCollision(i);
		}
	}
	
	//Integrates a pressure force on all tetrahedron surfaces
	//TODO - Pressure force flips due to inversion of the Normal direction, there needs to be a catch for this.
	private void SolveExternalVoxelPressureForce(float dt, float pressure)
	{
		for (int i = 0; i < numTets; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				int idThis = tetIds[4 * i + j];
				// The 3 opposite vertices ids
				int id0 = tetIds[4 * i + TetrahedronData.volIdOrder[j][0]];
				int id1 = tetIds[4 * i + TetrahedronData.volIdOrder[j][1]];
				int id2 = tetIds[4 * i + TetrahedronData.volIdOrder[j][2]];

				Vector3 id1_minus_id0 = pos[id1] - pos[id0];
				Vector3 id2_minus_id0 = pos[id2] - pos[id0];
				Vector3 cross = Vector3.Cross(id1_minus_id0, id2_minus_id0);

				float faceArea = cross.magnitude * 0.5f;
				Vector3 normal = cross.normalized;

				float pressureForce = pressure * faceArea;

				// Apply pressure force to each vertex of the face
				if (invMass[id0] != 0)
				{
					vel[id0] += (pressureForce * invMass[id0]) * normal * dt;
				}
				if (invMass[id1] != 0)
				{
					vel[id1] += (pressureForce * invMass[id1]) * normal * dt;
				}
				if (invMass[id2] != 0)
				{
					vel[id2] += (pressureForce * invMass[id2]) * normal * dt;
				}
			}
		}
	}

	private void SolvePressureForce(float dt, float pressure, int[] voxIDs, int[] face)
    {
        for (int i = 0; i < voxIDs.Length; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                // The id's of all particles on the face
				int[] vertexMapping = tetraData.GetVertexMapping;

				int id0 = vertexMapping[8 * voxIDs[i] + face[0]];
                int id1 = vertexMapping[8 * voxIDs[i] + face[1]];
                int id2 = vertexMapping[8 * voxIDs[i] + face[2]];
                int id3 = vertexMapping[8 * voxIDs[i] + face[3]];

				//Something about this does not work for all faces due to the way in which the triangles are formed.
                Vector3 id0_minus_id1 = pos[id0] - pos[id1];
                Vector3 id2_minus_id1 = pos[id2] - pos[id1];

				Vector3 id0_minus_id3 = pos[id0] - pos[id3];
				Vector3 id2_minus_id3 = pos[id2] - pos[id3];


                Vector3 crossF1 = Vector3.Cross(id0_minus_id1, id2_minus_id1);
                Vector3 crossF2 = Vector3.Cross(id0_minus_id3, id2_minus_id3);

                float faceAreaF1 = crossF1.magnitude * 0.5f;
                float faceAreaF2 = crossF2.magnitude * 0.5f;

                Vector3 normal = (crossF1.normalized-crossF2.normalized).normalized;

				/*Debug.DrawRay(pos[id0], -normal, Color.blue);
				Debug.DrawRay(pos[id1], -normal, Color.blue);
				Debug.DrawRay(pos[id2], -normal, Color.blue);
				Debug.DrawRay(pos[id3], -normal, Color.blue);*/

                float pressureForce = (pressure * (faceAreaF1+faceAreaF2))/4f;

                // Apply pressure force to each vertex of the face
                if (invMass[id0] != 0)
                {
                    vel[id0] += (pressureForce * invMass[id0]) * normal * dt;
                }
                if (invMass[id1] != 0)
                {
                    vel[id1] += (pressureForce * invMass[id1]) * normal * dt;
                }
                if (invMass[id2] != 0)
                {
                    vel[id2] += (pressureForce * invMass[id2]) * normal * dt;
                }
				if (invMass[id2] != 0)
                {
                    vel[id3] += (pressureForce * invMass[id2]) * normal * dt;
                }
            }
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

	
	private void debugLog()
	{
		//To calculate simulated displacement.
		/*Debug.Log("disps = " + (pos[beamLowerDisplacementPoss[0]].y- startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[1]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[2]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[3]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[4]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[5]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[6]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[7]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[8]].y - startingVerticalDisplacement));*/

		/*Debug.Log("disps = " + (pos[beamLowerDisplacementPoss[0]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[1]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[2]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[3]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[4]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[5]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[6]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[7]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[8]].x));*/

		//To calculate difference between real and simulated beam.
		/*Debug.Log("disps = " + (pos[beamLowerDisplacementPoss[0]].y- startingVerticalDisplacement - beamLowerDisplacementReal[0])
		+ " | " + (pos[beamLowerDisplacementPoss[1]].y - startingVerticalDisplacement - beamLowerDisplacementReal[1])
		+ " | " + (pos[beamLowerDisplacementPoss[2]].y - startingVerticalDisplacement - beamLowerDisplacementReal[2])
		+ " | " + (pos[beamLowerDisplacementPoss[3]].y - startingVerticalDisplacement - beamLowerDisplacementReal[3])
		+ " | " + (pos[beamLowerDisplacementPoss[4]].y - startingVerticalDisplacement - beamLowerDisplacementReal[4])
		+ " | " + (pos[beamLowerDisplacementPoss[5]].y - startingVerticalDisplacement - beamLowerDisplacementReal[5])
		+ " | " + (pos[beamLowerDisplacementPoss[6]].y - startingVerticalDisplacement - beamLowerDisplacementReal[6])
		+ " | " + (pos[beamLowerDisplacementPoss[7]].y - startingVerticalDisplacement - beamLowerDisplacementReal[7])
		+ " | " + (pos[beamLowerDisplacementPoss[8]].y - startingVerticalDisplacement - beamLowerDisplacementReal[8]));*/


		Debug.DrawRay(pos[cubeCompressionVoxels[0]*8], gravity, Color.yellow);
		Debug.DrawRay(pos[cubeCompressionVoxels[1]*8], gravity, Color.green);
		Debug.DrawRay(pos[cubeCompressionVoxels[2]*8], gravity, Color.blue);
		Debug.DrawRay(pos[cubeCompressionVoxels[3]*8], gravity, Color.red);
		Debug.DrawRay(pos[cubeCompressionVoxels[4]*8], gravity, Color.cyan);
		Debug.DrawRay(pos[cubeCompressionVoxels[23]*8], gravity, Color.gray);


		


		/*Debug.DrawRay(pos[beamLowerDisplacementPoss[0]], gravity, Color.yellow);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[1]], gravity, Color.green);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[2]], gravity, Color.red);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[3]], gravity, Color.blue);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[4]], gravity, Color.blue);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[5]], gravity, Color.gray);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[6]], gravity, Color.blue);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[7]], gravity, Color.cyan);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[8]], gravity, Color.red);*/

		/*Debug.DrawRay(pos[beamLowerDisplacementPoss[0]], gravity, Color.yellow);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[1]], gravity, Color.green);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[2]], gravity, Color.red);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[3]], gravity, Color.blue);
		Debug.DrawRay(pos[beamLowerDisplacementPoss[4]], gravity, Color.cyan);*/

		/*Debug.Log("disps = " + (pos[beamLowerDisplacementPoss[0]].y- startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[1]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[2]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[3]].y - startingVerticalDisplacement)
		+ " | " + (pos[beamLowerDisplacementPoss[4]].y - startingVerticalDisplacement));*/

		/*Debug.Log("disps = " + (pos[cube[0]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[1]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[2]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[3]].x)
		+ " | " + (pos[beamLowerDisplacementPoss[4]].x));*/
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

