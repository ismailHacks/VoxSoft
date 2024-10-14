using System.Collections;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.UI;


//Same as SoftBodySimulation but is using Vector3 instead of arrays where an index in the array is x, y, or z 
public class SoftBodySimulationVectors : IGrabbable
{
	//Tetrahedralizer data structures
	private readonly TetrahedronData tetraData;
	public voxelTet vox;
	private readonly int[] tetIds;
	private readonly int[] tetEdgeIds;
	private readonly float density = 1000; //kg/m^3
	private float startingVerticalDisplacement = 0.15f;
	public bool converged = false;

	/*public static int[] upperForce = new int[] {0,1,2,3,4,5,10,11,12,13,14,15,20,21,22,23,24,25,30,31,32,33,34,35,40,41,42,43,44,45,50,51,52,53};
	public static int[] rightForce = new int[] {2,3,4,5,12,13,14,15,16,17,18,19,21,22,25,26,34,35,36,37,44,45,46,47,48,49,50,51,53,54,57,58,66,67,68,69,76,77,78,80,81,82,85,86,89,90};
	public static int[] leftForce = new int[rightForce.Length];
	public static int[] lowerForce = new int[] {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
	public static int[] backForce = new int[] {4,5,6,7,24,25,26,27,44,45,46,47};
	public static int[] frontForce = new int[] {14,15,16,17,34,35,36,37,56,57,58,59,60,61,62,63};
	public static int[] lockFaceFront = new int[] {0,2,4,6,96,98,100,102,192,194,258,275};*/

	public static int[] forcePositiveX = new int[] {29,30,31,33,34,35,37,38,39};
	public static int[] forceNegativeX = new int[] {68,69,70,71,72,73,74,75,76};

	public static int[] forcePositiveY = new int[] {89,90,91,92,93,94,95,96,97};
	public static int[] forceNegativeY = new int[] {80,81,82,83,84,85,86,87,88};

	public static int[] forcePositiveZ = new int[] {50,51,52,54,55,56,58,59,60};
	public static int[] forceNegativeZ = new int[] {6,7,8,11,12,13,16,17,18};

	public static int[] lockFaceFront = new int[] {0,1,2,3,4,25,26,27,28,45,46,47,48,65,66,67,80,81,82,83,84,85,86,87,88};
	public static int[] lockFaceFront2 = new int[] {0};




	private readonly Vector3[] pos;
	private readonly Vector3[] prevPos;
	private readonly Vector3[] vel;

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

		/*for (int i = 0; i < rightForce.Length; i++)
		{
			leftForce[i] = rightForce[i]+96;
		}

		for (int i = 0; i < lowerForce.Length; i++)
		{
			lowerForce[i] = lowerForce[i]+258;
		}

		for (int i = 0; i < upperForce.Length; i++)
		{
			upperForce[i] = upperForce[i]+192;
		}

		for (int i = 0; i < backForce.Length; i++)
		{
			backForce[i] = backForce[i]+192;
		}

		for (int i = 0; i < frontForce.Length; i++)
		{
			frontForce[i] = frontForce[i]+192;
		}*/
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
		}

		float pMass = totalMass/numParticles;
		Debug.Log("Num Vertices = " + numParticles);
		Debug.Log("Mass = " + totalMass);

		for (int i = 0; i < numParticles; i++)
		{
			invMass[i] = 1/pMass;
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
		Simulate(dt, numSubSteps, edgeCompliance, volCompliance, dampingCoefficient, pressure);
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
		/*SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelPositiveX,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelPositiveX2,  Color.red);
		SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelNegativeX,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelNegativeX2,  Color.red);

		SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelPositiveY,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelPositiveY2,  Color.red);
		SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelNegativeY,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelNegativeY2,  Color.red);

		SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelPositiveZ,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelPositiveZ2,  Color.red);
		SolvePressureForce(dt, pressure, forcePositiveX2, voxelTet.voxelNegativeZ,  Color.green);
		SolvePressureForce(dt, pressure, forcePositiveX3, voxelTet.voxelNegativeZ2,  Color.red);*/

		SolvePressureForce(dt, pressure, forceNegativeY, voxelTet.voxelNegativeY);
		SolvePressureForce(dt, pressure, forcePositiveY, voxelTet.voxelPositiveY);
		SolvePressureForce(dt, pressure, forceNegativeZ, voxelTet.voxelNegativeZ);
		SolvePressureForce(dt, pressure, forcePositiveZ, voxelTet.voxelPositiveZ);
		SolvePressureForce(dt, pressure, forceNegativeX, voxelTet.voxelNegativeX);
		SolvePressureForce(dt, pressure, forcePositiveX, voxelTet.voxelPositiveX);

		/*SolvePressureForce2(dt, pressure, forceNegativeY, voxelTet.voxelNegativeY2);
		SolvePressureForce2(dt, pressure, forcePositiveY, voxelTet.voxelPositiveY2);
		SolvePressureForce2(dt, pressure, forceNegativeZ, voxelTet.voxelNegativeZ2);
		SolvePressureForce2(dt, pressure, forcePositiveZ, voxelTet.voxelPositiveZ2);
		SolvePressureForce2(dt, pressure, forceNegativeX, voxelTet.voxelNegativeX2);
		SolvePressureForce2(dt, pressure, forcePositiveX, voxelTet.voxelPositiveX2);*/
		

		lockFaces(lockFaceFront, voxelTet.voxelPositiveY);
		lockFaces2(lockFaceFront, voxelTet.voxelPositiveY2);


		forceMove(dt, dampingCoefficient);
		SolveEdges(dt, edgeCompliance);
		SolveVolumes(dt, volCompliance);
	}

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

	//TODO: This method is the bottleneck
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

	private void lockFaces2(int[] voxIDs, int[] face)
	{
        for (int i = 0; i < voxIDs.Length; i++)
        {
			for (int j = 0; j < 4; j++)
            {
				int[] vertexMapping = tetraData.GetVertexMapping;
                // The id's of all particles on the face
                invMass[vertexMapping[8 * (voxIDs[i]+98) + face[0]]] = 0f;
                invMass[vertexMapping[8 * (voxIDs[i]+98) + face[1]]] = 0f;
                invMass[vertexMapping[8 * (voxIDs[i]+98) + face[2]]] = 0f;
                invMass[vertexMapping[8 * (voxIDs[i]+98) + face[3]]] = 0f;
			}
		}
		int[] vertexMapping2 = tetraData.GetVertexMapping;
		Debug.DrawRay(pos[vertexMapping2[8 * (voxIDs[0]+98) + face[0]]], gravity, Color.black);
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
                int id1 = vertexMapping[8 * voxIDs[i] + face[1]];//This is free
                int id2 = vertexMapping[8 * voxIDs[i] + face[2]];
                int id3 = vertexMapping[8 * voxIDs[i] + face[3]];//This is free

				//Something about this does not work for all faces due to the way in which the triangles are formed.
                Vector3 id0_minus_id1 = pos[id0] - pos[id1];
                Vector3 id2_minus_id1 = pos[id2] - pos[id1];

				Vector3 id0_minus_id3 = pos[id0] - pos[id3];
				Vector3 id2_minus_id3 = pos[id2] - pos[id3];


                Vector3 crossF1 = Vector3.Cross(id0_minus_id1, id2_minus_id1);
                Vector3 crossF2 = Vector3.Cross(id0_minus_id3, id2_minus_id3);

                float faceAreaF1 = crossF1.magnitude * 0.5f;
                float faceAreaF2 = crossF2.magnitude * 0.5f;

                //Vector3 normal = (crossF1.normalized-crossF2.normalized).normalized;
                Vector3 normalF1 = crossF1.normalized;
                Vector3 normalF2 = crossF1.normalized;

				/*Debug.DrawRay(pos[id0], normalF1, color);
				Debug.DrawRay(pos[id1], normalF1, color);
				Debug.DrawRay(pos[id2], normalF1, color);
				Debug.DrawRay(pos[id3], normalF1, color);*/

                float pressureForceF1 = (pressure * faceAreaF1)/3f;
                float pressureForceF2 = (pressure * faceAreaF2)/3f;


                // Apply pressure force to each vertex of the face
                if (invMass[id0] != 0)
                {
                    vel[id0] += (pressureForceF1 * invMass[id0]) * normalF1 * dt;
                }
                if (invMass[id1] != 0)
                {
                    vel[id1] += (pressureForceF1 * invMass[id1]) * normalF1 * dt;
                }
                if (invMass[id2] != 0)
                {
                    vel[id2] += (pressureForceF1 * invMass[id2]) * normalF1 * dt;
                }

				// Apply pressure force to each vertex of the face
                if (invMass[id0] != 0)
                {
                    vel[id0] += (pressureForceF2 * invMass[id0]) * normalF2 * dt;
                }
                if (invMass[id2] != 0)
                {
                    vel[id2] += (pressureForceF2 * invMass[id1]) * normalF2 * dt;
                }
                if (invMass[id3] != 0)
                {
                    vel[id3] += (pressureForceF2 * invMass[id2]) * normalF2 * dt;
                }

            }
        }
    }

	private void SolvePressureForce2(float dt, float pressure, int[] voxIDs, int[] face)
    {
        for (int i = 0; i < voxIDs.Length; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                // The id's of all particles on the face
				int[] vertexMapping = tetraData.GetVertexMapping;

				int id0 = vertexMapping[8 * (voxIDs[i]+98) + face[0]];
                int id1 = vertexMapping[8 * (voxIDs[i]+98) + face[1]];//This is free
                int id2 = vertexMapping[8 * (voxIDs[i]+98) + face[2]];
                int id3 = vertexMapping[8 * (voxIDs[i]+98) + face[3]];//This is free

				//Something about this does not work for all faces due to the way in which the triangles are formed.
                Vector3 id0_minus_id1 = pos[id0] - pos[id1];
                Vector3 id2_minus_id1 = pos[id2] - pos[id1];

				Vector3 id0_minus_id3 = pos[id0] - pos[id3];
				Vector3 id2_minus_id3 = pos[id2] - pos[id3];


                Vector3 crossF1 = Vector3.Cross(id0_minus_id1, id2_minus_id1);
                Vector3 crossF2 = Vector3.Cross(id0_minus_id3, id2_minus_id3);

                float faceAreaF1 = crossF1.magnitude * 0.5f;
                float faceAreaF2 = crossF2.magnitude * 0.5f;

                //Vector3 normal = (crossF1.normalized-crossF2.normalized).normalized;
                Vector3 normalF1 = crossF1.normalized;
                Vector3 normalF2 = crossF1.normalized;

				/*Debug.DrawRay(pos[id0], normalF1, color);
				Debug.DrawRay(pos[id1], normalF1, color);
				Debug.DrawRay(pos[id2], normalF1, color);
				Debug.DrawRay(pos[id3], normalF1, color);*/

                float pressureForceF1 = (pressure * faceAreaF1)/3f;
                float pressureForceF2 = (pressure * faceAreaF2)/3f;


                // Apply pressure force to each vertex of the face
                if (invMass[id0] != 0)
                {
                    vel[id0] += (pressureForceF1 * invMass[id0]) * normalF1 * dt;
                }
                if (invMass[id1] != 0)
                {
                    vel[id1] += (pressureForceF1 * invMass[id1]) * normalF1 * dt;
                }
                if (invMass[id2] != 0)
                {
                    vel[id2] += (pressureForceF1 * invMass[id2]) * normalF1 * dt;
                }

				// Apply pressure force to each vertex of the face
                if (invMass[id0] != 0)
                {
                    vel[id0] += (pressureForceF2 * invMass[id0]) * normalF2 * dt;
                }
                if (invMass[id2] != 0)
                {
                    vel[id2] += (pressureForceF2 * invMass[id1]) * normalF2 * dt;
                }
                if (invMass[id3] != 0)
                {
                    vel[id3] += (pressureForceF2 * invMass[id2]) * normalF2 * dt;
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
		int[] vertexMapping = tetraData.GetVertexMapping;

		//4 is the index position on the specific voxel. *8 is for moving across the voxel id and then the number after is the voxelID to access
		/*Debug.DrawRay(pos[vertexMapping[4 + 8 * 89]], gravity, Color.blue);
		Debug.DrawRay(pos[vertexMapping[4 + 8 * 90]], gravity, Color.red);
		Debug.DrawRay(pos[vertexMapping[4 + 8 * 91]], gravity, Color.green);
		Debug.DrawRay(pos[vertexMapping[4 + 8 * 92]], gravity, Color.yellow);*/
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

