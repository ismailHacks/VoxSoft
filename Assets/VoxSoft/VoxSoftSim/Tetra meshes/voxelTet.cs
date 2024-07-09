//Using this test to add tetrahedrons in a voxel like way
//Uses 5 tetrahedrons per voxel

using System.Collections;
using System.Collections.Generic;
using JetBrains.Annotations;
using Unity.Mathematics;
using UnityEngine;
using Unity.Jobs;
using System.Security.Principal;
using System.Threading.Tasks;
using Unity.Collections;
using Unity.Burst;

public class voxelTet : TetrahedronData
{
	//Have to make sure number of voxels is correct to what is actually created!
	private static int noVoxels = 292;
	public static float voxelScale;
	private int globalVoxelCount = 0;
	private int connectionCount = 0;
	
	public int[] vertexMapping = new int[8*noVoxels];
	public float[] vertsVoxelMesh = new float[24*noVoxels];
	public int[] tetIdsVoxelMesh = new int[20*noVoxels];
	private int[] tetEdgeIdsVoxelMesh = new int[36*noVoxels];
	private int[] tetSurfaceTriIdsVoxelMesh = new int[48*noVoxels];

	//Getters
	public override float[] GetVerts => vertsVoxelMesh;
	public override int[] GetTetIds => tetIdsVoxelMesh;
	public override int[] GetTetEdgeIds => tetEdgeIdsVoxelMesh;
	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIdsVoxelMesh;
	public override int[] GetVertexMapping => vertexMapping;

	public voxelTet(float scale)
	{
		voxelScale = scale;
		float startTime = Time.realtimeSinceStartup;

		//Right Side
		makeCuboid(0,30,0,2,4,1);
		makeCuboid(2,30,0,4,6,1);
		makeCuboid(6,30,0,2,4,1);
		makeCuboid(8,30,0,4,6,1);
		makeCuboid(12,30,0,2,4,1);
		makeCuboid(14,30,0,4,6,1);
		//Left Side
		makeCuboid(0,30,3,2,4,1);
		makeCuboid(2,30,3,4,6,1);
		makeCuboid(6,30,3,2,4,1);
		makeCuboid(8,30,3,4,6,1);
		makeCuboid(12,30,3,2,4,1);
		makeCuboid(14,30,3,4,6,1);
		//Upper First
		makeCuboid(0,33,1,2,1,2);
		makeCuboid(2,33,1,1,3,2);
		makeCuboid(3,35,1,2,1,2);
		makeCuboid(5,33,1,1,3,2);
		//Upper Second
		makeCuboid(6,33,1,2,1,2);
		makeCuboid(8,33,1,1,3,2);
		makeCuboid(9,35,1,2,1,2);
		makeCuboid(11,33,1,1,3,2);
		//Upper Third
		makeCuboid(12,33,1,2,1,2);
		makeCuboid(14,33,1,1,3,2);
		makeCuboid(15,35,1,2,1,2);
		makeCuboid(17,30,1,1,6,2);
		//Lower
		makeCuboid(0,30,1,17,1,2);

		Debug.Log("Number of Voxels = " + globalVoxelCount);
		//Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
		combineVoxels(startTime);
		//Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
	}

	//
	// Design Library
	//

	private void makeActuator(int posX, int posY, int posZ, 
	float width, float wallThickness, float capHeight, 
	float totalHeight)
	{
		makeCylinder(posX,posY,posZ,width,capHeight);
		makeTube(posX,(int)capHeight,posZ,width,width-wallThickness,totalHeight-2*capHeight);
		makeCylinder(posX,(int)capHeight+(int)(totalHeight-2*capHeight),posZ,width,capHeight);
	}

	//
	// Core Shape Library
	//

	private void makeGyroid(int posX, int posY, int posZ, float width, float sensitivity)
	{
		for (int i = 0; i < (int)width; i++)
		{
			for (int j = (int)-0; j < width; j++)
			{
				for (int k = (int)-width/2; k < width/2; k++)
				{
					float gyroid = Mathf.Sin(i)*Mathf.Cos(j) + Mathf.Sin(j)*Mathf.Cos(k) + Mathf.Sin(k)*Mathf.Cos(i);
					if(gyroid < sensitivity && gyroid > -sensitivity)
					{
						makeVoxel(i+posX,j+posY,k+posZ);
					}
				}
			}	
		}
	}

	private void makeCylinder(int posX, int posY, int posZ, float radius, float height)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = (int)-radius; i < radius; i++)
			{
				for (int j = (int)-radius; j < radius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<radius)
					{
						makeVoxel(i+posX,k+posY,j+posZ);
					}
				}
			}	
		}
	}

	private void makeTube(int posX, int posY, int posZ, float outerRadius, float innerRadius, float height)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = (int)-outerRadius; i < outerRadius; i++)
			{
				for (int j = (int)-outerRadius; j < outerRadius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<outerRadius && Mathf.Sqrt(i*i+j*j)>innerRadius)
					{
						makeVoxel(i+posX,k+posY,j+posZ);
					}
				}
			}	
		}
	}

	private void makeCuboid(int posX, int posY, int posZ, int length, int height, int width)
	{
		for (int k = 0; k < height; k++)
		{
			for (int j = 0; j < width; j++)
			{
				for (int i = 0; i < length; i++)
				{
					makeVoxel(posX+i, posY+k, posZ+j);
				}
			}
		}
	}
	
	//
	// Voxel Generation
	//

	private void makeVoxel(int posX, int posY, int posZ)
	{
		for (int i = 0; i < verts.Length/3; i++)
		{
			vertsVoxelMesh[3*i+(24*globalVoxelCount)] = (verts[3*i]+posX)*voxelScale;
			vertsVoxelMesh[3*i+1+(24*globalVoxelCount)] = (verts[3*i+1]+posY)*voxelScale;
			vertsVoxelMesh[3*i+2+(24*globalVoxelCount)] = (verts[3*i+2]+posZ)*voxelScale;
		}

		for (int i = 0; i < tetIds.Length; i++)
		{
			tetIdsVoxelMesh[i+20*globalVoxelCount] = tetIds[i]+8*globalVoxelCount;
		}

		for (int i = 0; i < tetEdgeIds.Length; i++)
		{
			tetEdgeIdsVoxelMesh[i+36*globalVoxelCount] = tetEdgeIds[i]+8*globalVoxelCount;
		}

		for (int i = 0; i < tetSurfaceTriIds.Length; i++)
		{
			tetSurfaceTriIdsVoxelMesh[i+48*globalVoxelCount] = tetSurfaceTriIds[i]+8*globalVoxelCount;
		}
		globalVoxelCount++;
	}

	private void combineVoxels(float startTime)
	{
		HashSet<int> processedIndices = new HashSet<int>();
		int numVertices = vertsVoxelMesh.Length / 3;
		vertexMapping = new int[numVertices];

		for (int v = 0; v < numVertices; v++)
		{
			vertexMapping[v] = v;
		}

		//Debug.Log("VB = " + vertexMapping[2]);

		for (int i = 0; i < vertsVoxelMesh.Length / 3; i++)
		{
			if (processedIndices.Contains(i))
			{
				continue; // Skip already processed indices
			}
			int internalCount = 0;
			Vector3 iPosition = new Vector3(vertsVoxelMesh[3 * i], vertsVoxelMesh[3 * i + 1], vertsVoxelMesh[3 * i + 2]);

			for (int j = i + 1; j < vertsVoxelMesh.Length / 3; j++)
			{
				Vector3 jPosition = new Vector3(vertsVoxelMesh[3 * j], vertsVoxelMesh[3 * j + 1], vertsVoxelMesh[3 * j + 2]);
				if (iPosition == jPosition)
				{
					processedIndices.Add(j);
					connectionCount++;
					internalCount++;

					if(internalCount > 7)
					{
						internalCount = 0;
						continue;
					}
					
					for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
					{
						if (tetEdgeIdsVoxelMesh[k] == j)
						{
							tetEdgeIdsVoxelMesh[k] = i;
						}
					}

					for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
					{
						if (tetIdsVoxelMesh[l] == j)
						{
							tetIdsVoxelMesh[l] = i;
						}
					}

					for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
					{
						if (tetSurfaceTriIdsVoxelMesh[m] == j)
						{
							tetSurfaceTriIdsVoxelMesh[m] = i;
						}
					}
					vertexMapping[j] = i;
				}
			}
		}
		//Debug.Log("VA = " + vertexMapping[11]);
		//Debug.Log("CC = " + connectionCount + " - GVC = " + globalVoxelCount);
	}

	//Vertices (x, y, z) for a tetrahedral voxel
	//Provides the vertices of each particle in the voxel
	private float[] verts =
	{
		0,0,1,
		1,1,1,
		1,0,0,
		0,1,0,
		0,0,0,
		1,1,0,
		1,0,1,
		0,1,1
	};

	//Provides the ID position of the vertices that make up the tetrahedral voxel
	private int[] tetIds =
	{
		0,1,2,3,
		0,3,2,4,
		1,2,3,5,
		0,2,1,6,
		0,1,3,7
	};

	//Provides the connections between each one of the edges in the tetrahedral voxel
	//unlike tetIds the edges should not be repeated with connecting tetrahedrals as they will be looped over.
	private int[] tetEdgeIds =
	{
		//0,2, 2,1, 1,0, 1,3, 3,0, 2,3
		0,1, 1,2, 2,0, 0,3, 1,3, 2,3, //Centre Tetrahedron
		0,4, 2,4, 3,4, //Outer Tetrahedron 1
		1,5, 2,5, 3,5, //Outer Tetrahedron 2
		0,6, 1,6, 2,6,  //Outer Tetrahedron 3
		0,7, 1,7, 3,7  //Outer Tetrahedron 4
	};

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		//0,1,3, 1,2,3, 0,3,2, 0,2,1, //Inner Tetrahedron - Don't need to render inner tetrahedron as it is covered by the outer tetrahedrons
		0,3,4, 4,3,2, 0,4,2, 2,3,0, //Outer Tetrahedron 1
		2,3,5, 2,5,1, 3,1,5, 1,3,2, //Outer Tetrahedron 2
		0,6,1, 0,2,6, 1,6,2, 0,1,2, //Outer Tetrahedron 3
		0,7,3, 0,1,7, 1,3,7, 0,3,1  //Outer Tetrahedron 4
	};

    public static int[] voxelRight = new int[] {1, 6, 2, 5};
    public static int[] voxelLeft = new int[] {0, 7, 3, 4};
    public static int[] voxelFront= new int[] {1, 7, 0, 6};
    public static int[] voxelBack= new int[] {3, 5, 2, 4};
    public static int[] voxelTop= new int[] {1, 5, 3, 7};
    public static int[] voxelBottom= new int[] {0, 4, 2, 6};

	public int[] vertexMap()
	{
		return vertexMapping;
	}

	private void combineVoxelsLegacy() //Legacy code that has been superseeded by combineAndOptimiseVoxels()
	{
		for (int i = 0; i < vertsVoxelMesh.Length/3; i++)
		{
			for (int j = 0; j < vertsVoxelMesh.Length/3; j++)
			{
				if(vertsVoxelMesh[3*i]==vertsVoxelMesh[3*j] && vertsVoxelMesh[3*i+1]==vertsVoxelMesh[3*j+1] && vertsVoxelMesh[3*i+2]==vertsVoxelMesh[3*j+2] && i!=j)
				{
					connectionCount++;
					
					for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
					{
						if (tetEdgeIdsVoxelMesh[k] == j)
						{
							tetEdgeIdsVoxelMesh[k] = i;
						}
					}

					for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
					{
						if (tetIdsVoxelMesh[l] == j)
						{
							tetIdsVoxelMesh[l] = i;
						}
					}

					for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
					{
						if (tetSurfaceTriIdsVoxelMesh[m] == j)
						{
							tetSurfaceTriIdsVoxelMesh[m] = i;
						}
					}
				}
			}
		}
		Debug.Log("CC = " + connectionCount);
	}
}
