using System.Collections.Generic;
using UnityEngine;

public class voxelTet : TetrahedronData
{
	//Have to make sure number of voxels is correct to what is actually created!
	private static int cubicSize = 10;
	private static int noVoxels = 488;
	public static float voxelScale;
	private int globalVoxelCount = 0;

	public int[] vertexMapping = new int[8*noVoxels];
	public float[] vertsVoxelMesh = new float[24*noVoxels];
	public int[] tetIdsVoxelMesh = new int[20*noVoxels];
	private int[] tetEdgeIdsVoxelMesh = new int[36*noVoxels];
	private int[] tetSurfaceTriIdsVoxelMesh = new int[48*noVoxels];

	public bool[,,] voxelData = new bool[cubicSize, cubicSize, cubicSize];
	public Dictionary<Vector3Int, int> voxelPositionToID = new Dictionary<Vector3Int, int>();
	public Dictionary<string, List<int>> faceDirectionToVoxelIDs;
	public VoxelEnclosedSpaceDetector detector;

	//Getters
	public override float[] GetVerts => vertsVoxelMesh;
	public override int[] GetTetIds => tetIdsVoxelMesh;
	public override int[] GetTetEdgeIds => tetEdgeIdsVoxelMesh;
	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIdsVoxelMesh;
	public override int[] GetVertexMapping => vertexMapping;
	private int voxelID = 0;

	public voxelTet(float scale)
	{
		voxelScale = scale;
		float startTime = Time.realtimeSinceStartup;

		//makeCylinder(4,0,4,3,3,true);
		//GenerateRandomVoxelGrid(cubicSize);

		makeCuboid(0,0,0,10,10,10,true);
		makeCuboid(1,1,1,8,8,8,false);

		//makeCuboid(4,4,4,3,3,3,true);
		//makeCuboid(5,5,5,1,1,1,false);

		positionVoxels(voxelData);
		detector = new VoxelEnclosedSpaceDetector();
    	faceDirectionToVoxelIDs = detector.DetectEnclosedSpaces(voxelData, voxelPositionToID);

		//Debug.Log("Number of Voxels = " + globalVoxelCount);
		combineVoxels();
		//Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
	}

	private void positionVoxels(bool[,,] voxelData)
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
						globalVoxelCount++;
                    	voxelID++;
					}
				}
			}
		}
	}

	public void GenerateRandomVoxelGrid(int N)
	{
		System.Random random = new System.Random();

		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				for (int z = 0; z < N; z++)
				{
					bool isActive = random.NextDouble() >= 0.2; // Randomly assign true or false
					if (isActive)
					{
						voxelData[x, y, z] = true;
					}
					else
					{
						voxelData[x, y, z] = false;
					}
				}
			}
		}
	}

	//
	// Actuator Design Library
	//

	private void makeCylindricalActuator(int posX, int posY, int posZ, 
	float width, float wallThickness, float capHeight, 
	float totalHeight)
	{
		makeCylinder(posX,posY,posZ,width,capHeight, true);
		makeTube(posX,(int)capHeight,posZ,width,width-wallThickness,totalHeight-2*capHeight, true);
		makeCylinder(posX,(int)capHeight+(int)(totalHeight-2*capHeight),posZ,width,capHeight, true);
	}

	private void makePneuflexActuator(int posX, int posY, int posZ, 
	int width, int maxHeight, int innerHeight,int innerLength, int cellLength, 
	int numberCells, int wallThickness)
	{
		for (int i = 0; i < numberCells; i++)
		{
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ,innerLength,innerHeight,wallThickness, true);
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ+width-wallThickness,innerLength,innerHeight,wallThickness, true);
			makeCuboid(posX+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,innerLength,wallThickness,width-wallThickness*2, true);
			makeCuboid(posX+i*(innerLength+cellLength),posY,posZ+wallThickness,innerLength,wallThickness,width-wallThickness*2, true);
			
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ,cellLength,maxHeight,wallThickness, true);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ+width-wallThickness,cellLength,maxHeight,wallThickness, true);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY+maxHeight-wallThickness,posZ+wallThickness,cellLength,wallThickness,width-wallThickness*2, true);
			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY,posZ+wallThickness,cellLength,wallThickness,width-wallThickness*2, true);

			makeCuboid(posX+innerLength+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,wallThickness,maxHeight-innerHeight,width-wallThickness*2, true);
			makeCuboid(posX+innerLength+cellLength-wallThickness+i*(innerLength+cellLength),posY+innerHeight-wallThickness,posZ+wallThickness,wallThickness,maxHeight-innerHeight,width-wallThickness*2, true);
		}
		makeCuboid(posX+numberCells*(innerLength+cellLength)-wallThickness,posY+wallThickness,posZ+wallThickness,wallThickness,innerHeight-wallThickness*2,width-wallThickness*2, true);
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

	private void makeGyroid(int posX, int posY, int posZ, float width, float sensitivity, float tp, bool addOrSubtract)
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
						voxelData[posX+i, posY+k, posZ+j] = addOrSubtract;
					}
				}
			}	
		}
	}

	private void makeCylinder(int posX, int posY, int posZ, float radius, float height, bool addOrSubtract)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = 0; i < 2*radius; i++)
			{
				for (int j = 0; j < 2*radius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<radius)
					{
						voxelData[posX+i, posY+k, posZ+j] = addOrSubtract;
					}
				}
			}	
		}
	}
	
	private void makeTube(int posX, int posY, int posZ, float outerRadius, float innerRadius, float height, bool addOrSubtract)
	{
		for (int k = 0; k < height; k++)
		{
			for (int i = (int)-outerRadius; i < outerRadius; i++)
			{
				for (int j = (int)-outerRadius; j < outerRadius; j++)
				{
					if(Mathf.Sqrt(i*i+j*j)<outerRadius && Mathf.Sqrt(i*i+j*j)>innerRadius)
					{
						voxelData[posX+i, posY+k, posZ+j] = addOrSubtract;
					}
				}
			}	
		}
	}

	private void makeCuboid(int posX, int posY, int posZ, int length, int height, int width, bool addOrSubtract)
	{
		for (int k = 0; k < height; k++)
		{
			for (int j = 0; j < width; j++)
			{
				for (int i = 0; i < length; i++)
				{
					voxelData[posX+i, posY+k, posZ+j] = addOrSubtract;
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
	// Voxel Generation
	//

	private void makeVoxel(int posX, int posY, int posZ)
	{
		for (int i = 0; i < verts.Length/3; i++)
		{
			//Debug.Log(vertsVoxelMesh.Length);
			//Debug.Log(vertsVoxelMesh[0]);
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
	}

	private void combineVoxels()
	{
		int numVertices = vertsVoxelMesh.Length / 3;
		vertexMapping = new int[numVertices];
		Dictionary<Vector3, int> positionToIndex = new Dictionary<Vector3, int>();
		Dictionary<Vector3, int> positionInternalCount = new Dictionary<Vector3, int>();

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

				internalCount++;
				positionInternalCount[iPosition] = internalCount;
				vertexMapping[i] = existingIndex;
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

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		0,3,4, 4,3,2, 0,4,2, 2,3,0, //Outer Tetrahedron 1
		2,3,5, 2,5,1, 3,1,5, 1,3,2, //Outer Tetrahedron 2
		0,6,1, 0,2,6, 1,6,2, 0,1,2, //Outer Tetrahedron 3
		0,7,3, 0,1,7, 1,3,7, 0,3,1  //Outer Tetrahedron 4
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
}
