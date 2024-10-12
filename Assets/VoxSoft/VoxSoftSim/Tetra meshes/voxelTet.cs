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
using System.Linq; // Add this to use the Sum() method

public class voxelTet : TetrahedronData
{
	//Have to make sure number of voxels is correct to what is actually created!
	private static int cubicSize = 7;
	//private static int noVoxels = cubicSize*cubicSize*cubicSize+cubicSize*2*cubicSize*2*cubicSize*2+cubicSize*5*cubicSize*5*cubicSize*5;
	private static int noVoxels = cubicSize*cubicSize*cubicSize;
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
	private int voxelID = 0;

	bool[,,] voxelData = new bool[cubicSize, cubicSize, cubicSize];
	Dictionary<Vector3Int, int> voxelPositionToID = new Dictionary<Vector3Int, int>();
	public Dictionary<string, List<int>> faceDirectionToVoxelIDs;

	public voxelTet(float scale)
	{
		voxelScale = scale;
		float startTime = Time.realtimeSinceStartup;

		/*for (int xPos = 0; xPos < cubicSize-3; xPos++)
		{
			for (int yPos = 0; yPos < cubicSize-3; yPos++)
			{
				for (int zPos = 0; zPos < cubicSize-3; zPos++)
				{
						voxelData[xPos, yPos, zPos] = true; // Filled
				}
			}
		}

		for (int xPos = 1; xPos < cubicSize - 4; xPos++)
		{
			for (int yPos = 1; yPos < cubicSize - 4; yPos++)
			{
				for (int zPos = 1; zPos < cubicSize - 4; zPos++)
				{
					voxelData[xPos, yPos, zPos] = false; // Empty inside
				}
			}
		}*/

		voxelData[0, 0, 0] = true;

		makeFreeActuator(voxelData);

		VoxelEnclosedSpaceDetector detector = new VoxelEnclosedSpaceDetector();
    	faceDirectionToVoxelIDs = detector.DetectEnclosedSpaces(voxelData, voxelPositionToID);

		/*foreach (var kvp in faceDirectionToVoxelIDs)
		{
			string faceDirection = kvp.Key;
			List<int> voxelIDs = kvp.Value;

			Debug.Log($"Face Direction: {faceDirection}");
			Debug.Log($"Voxel IDs: {string.Join(", ", voxelIDs)}");
		}*/

		Debug.Log("Number of Voxels = " + globalVoxelCount);
		combineVoxels();
		Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
	}

	//
	// Actuator Design Library
	//
	private void makeFreeActuator(bool[,,] voxelData)
	{
		int sizeX = voxelData.GetLength(0);
    	int sizeY = voxelData.GetLength(1);
    	int sizeZ = voxelData.GetLength(2);

		for (int xPos = 0; xPos < sizeX; xPos++)
		{
			for (int yPos = 0; yPos < sizeY; yPos++)
			{
				for (int zPos = 0; zPos < sizeZ; zPos++)
				{
					if (voxelData[xPos, yPos, zPos])
					{
						makeVoxel(xPos, yPos, zPos);
						Vector3Int position = new Vector3Int(xPos, yPos, zPos);
                    	voxelPositionToID[position] = voxelID;
                    	voxelID++;
					}
				}
			}
		}
	}



	private void makeCylindricalActuator(int posX, int posY, int posZ, 
	float width, float wallThickness, float capHeight, 
	float totalHeight)
	{
		makeCylinder(posX,posY,posZ,width,capHeight);
		makeTube(posX,(int)capHeight,posZ,width,width-wallThickness,totalHeight-2*capHeight);
		makeCylinder(posX,(int)capHeight+(int)(totalHeight-2*capHeight),posZ,width,capHeight);
	}

	private void makePneuflexActuator(int posX, int posY, int posZ, 
	int width, int maxHeight, int innerHeight,int innerLength, int cellLength, 
	int numberCells, int wallThickness)
	{
		for (int i = 0; i < numberCells; i++)
		{
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ,innerLength,innerHeight,wallThickness);
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ+width-wallThickness,innerLength,innerHeight,wallThickness);
			makeCuboid(posX+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,innerLength,wallThickness,width-wallThickness*2);
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ+wallThickness,innerLength,wallThickness,width-wallThickness*2);
			
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ,cellLength,maxHeight,wallThickness);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ+width-wallThickness,cellLength,maxHeight,wallThickness);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY+maxHeight-wallThickness,posZ+wallThickness,cellLength,wallThickness,width-wallThickness*2);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ+wallThickness,cellLength,wallThickness,width-wallThickness*2);

			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,wallThickness,maxHeight-innerHeight,width-wallThickness*2);
			makeCuboid(posX+innerLength+cellLength-wallThickness+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,wallThickness,maxHeight-innerHeight,width-wallThickness*2);
		}
		makeCuboid(posX+numberCells*(innerLength+cellLength)-wallThickness,posY+wallThickness,posZ+wallThickness,wallThickness,innerHeight-wallThickness*2,width-wallThickness*2);
	}

	private void makeSphericalActuator(int posX, int posY, int posZ, 
	float radius, float wallThickness)
	{
		for (int k = (int)-radius; k < radius; k++)
		{
			for (int j = (int)-radius; j < radius; j++)
			{
				for (int i = (int)-radius; i < radius; i++)
				{
					float distVal = Mathf.Sqrt(i*i+j*j+k*k);
					if (distVal < radius && distVal > (radius-wallThickness))
					{
						makeVoxel(posX+i, posY+k, posZ+j);
					}
				}
			}
		}
	}

	//
	// Core Shape Library
	//

	private void makeGyroid(int posX, int posY, int posZ, float width, float sensitivity, float tp)
	{
		for (int i = 0; i < (int)width; i++)
		{
			for (int j = (int)-0; j < width; j++)
			{
				for (int k = (int)-width/2; k < width/2; k++)
				{
					float gyroid = Mathf.Sin(i/tp)*Mathf.Cos(j/tp) + Mathf.Sin(j/tp)*Mathf.Cos(k/tp) + Mathf.Sin(k/tp)*Mathf.Cos(i/tp);
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
	
	//Currently gets messed up with environment collision
	private void makeSphere(int posX, int posY, int posZ, int radius)
	{
		for (int k = -radius; k < radius; k++)
		{
			for (int j = -radius; j < radius; j++)
			{
				for (int i = -radius; i < radius; i++)
				{
					float distVal = Mathf.Sqrt(i*i+j*j+k*k);
					if (distVal < radius)
					{
						makeVoxel(posX+i, posY+k, posZ+j);
					}
				}
			}
		}
	}

	//
	// Core Shape Library Alternate Mesh
	//
	private void makeGyroid2(int posX, int posY, int posZ, float width, float sensitivity, float tp)
	{
		for (int i = 0; i < (int)width; i++)
		{
			for (int j = (int)-0; j < width; j++)
			{
				for (int k = (int)-width/2; k < width/2; k++)
				{
					float gyroid = Mathf.Sin(i/tp)*Mathf.Cos(j/tp) + Mathf.Sin(j/tp)*Mathf.Cos(k/tp) + Mathf.Sin(k/tp)*Mathf.Cos(i/tp);
					if(gyroid < sensitivity && gyroid > -sensitivity)
					{
						makeVoxel2(i+posX,j+posY,k+posZ);
					}
				}
			}	
		}
	}

	private void makeCylinder2(int posX, int posY, int posZ, float radius, float height)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = (int)-radius; i < radius; i++)
			{
				for (int j = (int)-radius; j < radius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<radius)
					{
						makeVoxel2(i+posX,k+posY,j+posZ);
					}
				}
			}	
		}
	}

	private void makeTube2(int posX, int posY, int posZ, float outerRadius, float innerRadius, float height)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = (int)-outerRadius; i < outerRadius; i++)
			{
				for (int j = (int)-outerRadius; j < outerRadius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<outerRadius && Mathf.Sqrt(i*i+j*j)>innerRadius)
					{
						makeVoxel2(i+posX,k+posY,j+posZ);
					}
				}
			}	
		}
	}

	private void makeCuboid2(int posX, int posY, int posZ, int length, int height, int width)
	{
		for (int k = 0; k < height; k++)
		{
			for (int j = 0; j < width; j++)
			{
				for (int i = 0; i < length; i++)
				{
					makeVoxel2(posX+i, posY+k, posZ+j);
				}
			}
		}
	}

	//Currently gets messed up with environment collision
	private void makeSphere2(int posX, int posY, int posZ, int radius)
	{
		for (int k = -radius; k < radius; k++)
		{
			for (int j = -radius; j < radius; j++)
			{
				for (int i = -radius; i < radius; i++)
				{
					float distVal = Mathf.Sqrt(i*i+j*j+k*k);
					if (distVal < radius)
					{
						makeVoxel2(posX+i, posY+k, posZ+j);
					}
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

	private void makeVoxel2(int posX, int posY, int posZ)
	{
		for (int i = 0; i < verts.Length/3; i++)
		{
			vertsVoxelMesh[3*i+(24*globalVoxelCount)] = (verts[3*i]+posX)*voxelScale;
			vertsVoxelMesh[3*i+1+(24*globalVoxelCount)] = (verts[3*i+1]+posY)*voxelScale;
			vertsVoxelMesh[3*i+2+(24*globalVoxelCount)] = (verts[3*i+2]+posZ)*voxelScale;
		}

		for (int i = 0; i < tetIds.Length; i++)
		{
			tetIdsVoxelMesh[i+20*globalVoxelCount] = tetIds2[i]+8*globalVoxelCount;
		}

		for (int i = 0; i < tetEdgeIds.Length; i++)
		{
			tetEdgeIdsVoxelMesh[i+36*globalVoxelCount] = tetEdgeIds2[i]+8*globalVoxelCount;
		}

		for (int i = 0; i < tetSurfaceTriIds.Length; i++)
		{
			tetSurfaceTriIdsVoxelMesh[i+48*globalVoxelCount] = tetSurfaceTriIds2[i]+8*globalVoxelCount;
		}
		globalVoxelCount++;
	}

	private void combineVoxels()
	{
		int numVertices = vertsVoxelMesh.Length / 3;
		vertexMapping = new int[numVertices];
		Dictionary<Vector3, int> positionToIndex = new Dictionary<Vector3, int>();
		Dictionary<Vector3, int> positionInternalCount = new Dictionary<Vector3, int>();
		connectionCount = 0;

		for (int i = 0; i < numVertices; i++)
		{
			Vector3 iPosition = new Vector3(vertsVoxelMesh[3 * i], vertsVoxelMesh[3 * i + 1], vertsVoxelMesh[3 * i + 2]);

			if (positionToIndex.TryGetValue(iPosition, out int existingIndex))
			{
				// Check internalCount
				if (!positionInternalCount.TryGetValue(iPosition, out int internalCount))
				{
					internalCount = 0;
				}

				if (internalCount >= 7)
				{
					// Skip mapping this duplicate
					continue;
				}

				internalCount++;
				positionInternalCount[iPosition] = internalCount;

				vertexMapping[i] = existingIndex;
				connectionCount++;
			}
			else
			{
				positionToIndex[iPosition] = i;
				positionInternalCount[iPosition] = 0;
				vertexMapping[i] = i;
			}
		}

		// Now update the indices in the arrays
		for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
		{
			tetEdgeIdsVoxelMesh[k] = vertexMapping[tetEdgeIdsVoxelMesh[k]];
		}

		for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
		{
			tetIdsVoxelMesh[l] = vertexMapping[tetIdsVoxelMesh[l]];
		}

		for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
		{
			tetSurfaceTriIdsVoxelMesh[m] = vertexMapping[tetSurfaceTriIdsVoxelMesh[m]];
		}
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

	private int[] tetIds2 =
	{
		4,7,6,5,
		4,5,6,2,
		7,6,5,1,
		4,6,7,0,
		4,7,5,3
	};

	//Provides the connections between each one of the edges in the tetrahedral voxel
	//unlike tetIds the edges should not be repeated with connecting tetrahedrals as they will be looped over.
	private int[] tetEdgeIds =
	{
		0,1, 1,2, 2,0, 0,3, 1,3, 2,3, //Centre Tetrahedron
		0,4, 2,4, 3,4, //Outer Tetrahedron 1
		1,5, 2,5, 3,5, //Outer Tetrahedron 2
		0,6, 1,6, 2,6,  //Outer Tetrahedron 3
		0,7, 1,7, 3,7  //Outer Tetrahedron 4
	};

	private int[] tetEdgeIds2 =
	{
		4,7, 7,6, 6,4, 4,5, 7,5, 6,5,
		4,2, 6,2, 5,2,
		7,1, 6,1, 5,1,
		4,0, 7,0, 6,0,
		4,3, 7,3, 5,3
	};

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		0,3,4, 4,3,2, 0,4,2, 2,3,0, //Outer Tetrahedron 1
		2,3,5, 2,5,1, 3,1,5, 1,3,2, //Outer Tetrahedron 2
		0,6,1, 0,2,6, 1,6,2, 0,1,2, //Outer Tetrahedron 3
		0,7,3, 0,1,7, 1,3,7, 0,3,1  //Outer Tetrahedron 4
	};

	private int[] tetSurfaceTriIds2 =
	{
		4,5,2, 2,5,6, 4,2,6, 6,5,4,
		6,5,1, 6,1,7, 5,7,1, 7,5,6,
		4,0,7, 4,6,0, 7,0,6, 4,7,6,
		4,3,5, 4,7,3, 7,5,3, 4,5,7
	};

	public int[] vertexMap()
		{
			return vertexMapping;
		}
		
	//This gives the first face that is intersected when moving in the direction from outside the voxel.
    public static int[] voxelNegativeX = new int[] {1, 6, 2, 5};
    public static int[] voxelPositiveX = new int[] {0, 7, 3, 4};
    public static int[] voxelNegativeZ = new int[] {1, 7, 0, 6};
    public static int[] voxelPositiveZ = new int[] {3, 5, 2, 4};
    public static int[] voxelNegativeY = new int[] {1, 5, 3, 7};
    public static int[] voxelPositiveY = new int[] {0, 4, 2, 6};

	public static int[] voxelNegativeX2 = new int[] {7, 0, 6, 1};
    public static int[] voxelPositiveX2 = new int[] {4, 3, 5, 2};
    public static int[] voxelNegativeZ2 = new int[] {7, 3, 4, 0};
    public static int[] voxelPositiveZ2 = new int[] {5, 1, 6, 2};
    public static int[] voxelNegativeY2 = new int[] {7, 1, 5, 3};
    public static int[] voxelPositiveY2 = new int[] {4, 2, 6, 0};

    public static int[] voxelRight = new int[] {1, 6, 2, 5};
    public static int[] voxelLeft = new int[] {0, 7, 3, 4};
    public static int[] voxelFront = new int[] {1, 7, 0, 6};
    public static int[] voxelBack = new int[] {3, 5, 2, 4};
    public static int[] voxelTop = new int[] {1, 5, 3, 7};
    public static int[] voxelBottom = new int[] {0, 4, 2, 6};
}
