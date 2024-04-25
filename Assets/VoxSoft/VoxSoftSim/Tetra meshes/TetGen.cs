//Using this test to add tetrahedrons in a voxel like way
//Uses 5 tetrahedrons per voxel

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class tetGen : TetrahedronData
{
	private static int noTets = 4;
	private static float tetScale = 0.2f;
	private int globalTetCount = 0;


	public float[] vertsMesh = new float[12+3*(noTets-1)];
	public int[] tetIdsMesh = new int[4*noTets];
	private int[] tetEdgeIdsMesh = new int[12+6*(noTets-1)];
	private int[] tetSurfaceTriIdsMesh = new int[12+9*(noTets-1)];

	List<Vector2> ChangeID = new List<Vector2>();

	//Getters
	public override float[] GetVerts => vertsMesh;
	public override int[] GetTetIds => tetIdsMesh;
	public override int[] GetTetEdgeIds => tetEdgeIdsMesh;
	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIdsMesh;

	public tetGen()
	{
		float startTime = Time.realtimeSinceStartup;

		seedTet();
		makeTet(0,0,-1);
		makeTet(1,1,-1);
		makeTet(1,1,2);

		Debug.Log(globalTetCount);
		//Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
		//Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
	}

	private void seedTet()
	{
		for (int i = 0; i < verts.Length/3; i++)
		{
			vertsMesh[3*i] = verts[3*i]*tetScale;
			vertsMesh[3*i+1] = verts[3*i+1]*tetScale;
			vertsMesh[3*i+2] = verts[3*i+2]*tetScale;
		}

		for (int i = 0; i < tetIds.Length; i++)
		{
			tetIdsMesh[i] = tetIds[i];
		}

		for (int i = 0; i < tetEdgeIds.Length; i++)
		{
			tetEdgeIdsMesh[i] = tetEdgeIds[i];
		}

		for (int i = 0; i < tetSurfaceTriIds.Length; i++)
		{
			tetSurfaceTriIdsMesh[i] = tetSurfaceTriIds[i];
		}
	}

	private void makeTet(int posX, int posY, int posZ)
    {
		vertsMesh[12+3*globalTetCount] = posX*tetScale;
		vertsMesh[12+3*globalTetCount+1] = posY*tetScale;
		vertsMesh[12+3*globalTetCount+2] = posZ*tetScale;

        //Need to find the closest 3 vertices and make them the tetId points
        int max3ID, max2ID, max1ID;
        findClosestPlaneVertices(posX, posY, posZ, out max3ID, out max2ID, out max1ID);

        tetIdsMesh[4+4*globalTetCount] = 4+globalTetCount;
		tetIdsMesh[4+4*globalTetCount+1] = max3ID;
		tetIdsMesh[4+4*globalTetCount+2] = max2ID;
		tetIdsMesh[4+4*globalTetCount+3] = max1ID;

		tetEdgeIdsMesh[12+6*globalTetCount] = max3ID;
		tetEdgeIdsMesh[12+6*globalTetCount+1] = 4+globalTetCount;
		tetEdgeIdsMesh[12+6*globalTetCount+2] = max2ID;
		tetEdgeIdsMesh[12+6*globalTetCount+3] = 4+globalTetCount;
		tetEdgeIdsMesh[12+6*globalTetCount+4] = max1ID;
		tetEdgeIdsMesh[12+6*globalTetCount+5] = 4+globalTetCount;

		tetSurfaceTriIdsMesh[12+9*globalTetCount] = 4+globalTetCount;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+1] = max3ID;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+2] = max1ID;

		tetSurfaceTriIdsMesh[12+9*globalTetCount+3] = 4+globalTetCount;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+4] = max1ID;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+5] = max2ID;

		tetSurfaceTriIdsMesh[12+9*globalTetCount+6] = 4+globalTetCount;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+7] = max2ID;
		tetSurfaceTriIdsMesh[12+9*globalTetCount+8] = max3ID;
		globalTetCount++;
    }

    private void findClosestPlaneVertices(int posX, int posY, int posZ, out int max3ID, out int max2ID, out int max1ID)
    {
		float averageVertexDistance;
        max3ID = 0;
        max2ID = 0;
        max1ID = 0;
		float maxAverage = 0;
		int maxAverageID = 0;


        //Getting closest vertices
        //TODO Need to make sure it finds the vertises that are closest of a single face, rather than that from multiple faces

        for (int i = 0; i < ((12+9*globalTetCount)/3); i++)
        {
            Vector3 currentPos = new Vector3(posX * tetScale, posY * tetScale, posZ * tetScale);
			int vertexID1 = tetSurfaceTriIdsMesh[3*i];
			int vertexID2 = tetSurfaceTriIdsMesh[3*i+1];
			int vertexID3 = tetSurfaceTriIdsMesh[3*i+2];
            Vector3 vertexPos1 = new Vector3(vertsMesh[3*vertexID1], vertsMesh[3*vertexID1+1], vertsMesh[3*vertexID1+2]);
			Vector3 vertexPos2 = new Vector3(vertsMesh[3*vertexID2], vertsMesh[3*vertexID2+1], vertsMesh[3*vertexID2+2]);
			Vector3 vertexPos3 = new Vector3(vertsMesh[3*vertexID3], vertsMesh[3*vertexID3+1], vertsMesh[3*vertexID3+2]);
			averageVertexDistance = ((currentPos-vertexPos1).magnitude + (currentPos-vertexPos2).magnitude + (currentPos-vertexPos3).magnitude)/3;
			
			//float vertexDistance = Mathf.Abs((currentPos - measureToPlanePos).magnitude);
            Debug.Log("Av PlD" + i + " - " + averageVertexDistance);
			if (averageVertexDistance != 0)
			{
				if (averageVertexDistance < maxAverage || maxAverage == 0)
				{
					maxAverage = averageVertexDistance;
					maxAverageID = i;
				}
			}
			//Debug.Log("MaxAverage = " + maxAverage + " - ID = " + maxAverageID);
			max1ID = tetSurfaceTriIdsMesh[3*maxAverageID];
			max2ID = tetSurfaceTriIdsMesh[3*maxAverageID+1];
			max3ID = tetSurfaceTriIdsMesh[3*maxAverageID+2];
        }
		Debug.Log("Max3 = " + max3ID + " Max2 = " + max2ID + " Max1 = " + max1ID);
    }

    //Vertices (x, y, z) for the first tetrahedral
    private float[] verts =
	{
		0,0,1,
		0,1,0,
		1,0,0,
		0,0,0
	};

	//Provides the ID position of the vertices that make up the tetrahedral voxel
	private int[] tetIds =
	{
		0,1,2,3
	};

	//Provides the connections between each one of the edges in the tetrahedral voxel
	//unlike tetIds the edges should not be repeated with connecting tetrahedrals as they will be looped over.
	private int[] tetEdgeIds =
	{
		0,1, 1,2, 2,0, 0,3, 1,3, 2,3
	};

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		0,1,3, 1,2,3, 0,3,2, 0,2,1
	};
}
