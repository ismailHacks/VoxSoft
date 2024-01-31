//Using this test to add tetrahedrons in a voxel like way

using System.Collections;
using System.Collections.Generic;
using Unity.Mathematics;
using UnityEngine;
using Unity.Jobs;
using Unity.Collections;
using Unity.Burst;

public class voxelTetMultiThreaded : TetrahedronData
{
	private static int noVoxels = 1000;
	private static float voxelScale = 0.1f;
	private int globalVoxelCount = 0;
	private int connectionCount = 0;

	public float[] vertsVoxelMesh = new float[24*noVoxels];
	public int[] tetIdsVoxelMesh = new int[20*noVoxels];
	private int[] tetEdgeIdsVoxelMesh = new int[36*noVoxels];
	private int[] tetSurfaceTriIdsVoxelMesh = new int[48*noVoxels];

	//Getters
	public override float[] GetVerts => vertsVoxelMesh;
	public override int[] GetTetIds => tetIdsVoxelMesh;
	public override int[] GetTetEdgeIds => tetEdgeIdsVoxelMesh;
	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIdsVoxelMesh;

	public voxelTetMultiThreaded()
	{
		float startTime = Time.realtimeSinceStartup;
		makeActuator(0,0,0,4,2,3,8);
		Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
		//combineAndOptimizeVoxels();
		combineVoxelsMulti(); //Currently burst compiling
		Debug.Log(((Time.realtimeSinceStartup-startTime)*1000f)+" ms");
		Debug.Log(globalVoxelCount);
	}

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

	private void makeVoxelCombiner(int posX, int posY, int posZ)
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



		for (int i = 0; i < vertsVoxelMesh.Length/3; i++)
		{
			for (int j = 0; j < vertsVoxelMesh.Length/3; j++)
			{
				if(vertsVoxelMesh[3*i]==vertsVoxelMesh[3*j] && vertsVoxelMesh[3*i+1]==vertsVoxelMesh[3*j+1] && vertsVoxelMesh[3*i+2]==vertsVoxelMesh[3*j+2] && i!=j)
				{
					connectionCount++;
					int iPos = 3*i;
					int jPos = 3*j;
					
					for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
					{
						if (tetEdgeIdsVoxelMesh[k] == jPos/3)
						{
							tetEdgeIdsVoxelMesh[k] = iPos/3;
						}
					}

					for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
					{
						if (tetIdsVoxelMesh[l] == jPos/3)
						{
							tetIdsVoxelMesh[l] = iPos/3;
						}
					}

					for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
					{
						if (tetSurfaceTriIdsVoxelMesh[m] == jPos/3)
						{
							tetSurfaceTriIdsVoxelMesh[m] = iPos/3;
						}
					}
				}
			}
		}



	}

	private void combineAndOptimizeVoxels()
	{
		HashSet<int> processedIndices = new HashSet<int>();

		for (int i = 0; i < vertsVoxelMesh.Length / 3; i++)
		{
			if (processedIndices.Contains(i))
			{
				continue; // Skip already processed indices
			}

			Vector3 iPosition = new Vector3(vertsVoxelMesh[3 * i], vertsVoxelMesh[3 * i + 1], vertsVoxelMesh[3 * i + 2]);

			for (int j = i + 1; j < vertsVoxelMesh.Length / 3; j++)
			{
				Vector3 jPosition = new Vector3(vertsVoxelMesh[3 * j], vertsVoxelMesh[3 * j + 1], vertsVoxelMesh[3 * j + 2]);

				if (iPosition == jPosition)
				{

					processedIndices.Add(j);

					int iPos = 3 * i;
					int jPos = 3 * j;
					
					for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
					{
						if (tetEdgeIdsVoxelMesh[k] == jPos / 3)
						{
							tetEdgeIdsVoxelMesh[k] = iPos / 3;
						}
					}

					//Multithreading code
					/*NativeArray<int> tetIDSJobArrayExternal = new NativeArray<int>(tetIdsVoxelMesh.Length, Allocator.TempJob);

					for (int n = 0; n < tetIDSJobArrayExternal.Length; n++)
					{
						tetIDSJobArrayExternal[n] = tetIdsVoxelMesh[n];
					}

					tetIDSJob tetIDSJob = new tetIDSJob{tetIDSJobArrayInternal = tetIDSJobArrayExternal,i1=iPos,j1=jPos};
					JobHandle jobHandle = tetIDSJob.Schedule();
					jobHandle.Complete();

					for (int o = 0; o < tetIDSJobArrayExternal.Length; o++)
					{
						tetIdsVoxelMesh[o] = tetIDSJobArrayExternal[o];
					}
					tetIDSJobArrayExternal.Dispose();*/

					
					for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
					{
						if (tetIdsVoxelMesh[l] == jPos / 3)
						{
							tetIdsVoxelMesh[l] = iPos / 3;
						}
					}
					


					for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
					{
						if (tetSurfaceTriIdsVoxelMesh[m] == jPos / 3)
						{
							tetSurfaceTriIdsVoxelMesh[m] = iPos / 3;
						}
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

	private void makeActuator(int posX, int posY, int posZ, float width, float wallThickness, float capHeight, float totalHeight)
	{
		makeCylinder(posX,posY,posZ,width,capHeight);
		makeTube(posX,(int)capHeight,posZ,width,width-wallThickness,totalHeight-2*capHeight);
		makeCylinder(posX,(int)capHeight+(int)(totalHeight-2*capHeight),posZ,width,capHeight);
	}






	private void combineVoxelsMulti()
	{
		float startTime1= Time.realtimeSinceStartup;
		
		NativeArray<float> vertsJobArrayExternal = new NativeArray<float>(vertsVoxelMesh.Length, Allocator.TempJob);
		NativeArray<int> tetIDSJobArrayExternal = new NativeArray<int>(tetIdsVoxelMesh.Length, Allocator.TempJob);
		NativeArray<int> tetEdgeIDSJobArrayExternal = new NativeArray<int>(tetEdgeIdsVoxelMesh.Length, Allocator.TempJob);
		NativeArray<int> tetSurfaceTriIdsJobArrayExternal = new NativeArray<int>(tetSurfaceTriIdsVoxelMesh.Length, Allocator.TempJob);;
		//NativeHashSet<float> processedIndicesExternal;

		for (int n = 0; n < vertsJobArrayExternal.Length; n++)
		{
			vertsJobArrayExternal[n] = vertsVoxelMesh[n];
		}

		for (int n = 0; n < tetIDSJobArrayExternal.Length; n++)
		{
			tetIDSJobArrayExternal[n] = tetIdsVoxelMesh[n];
		}

		for (int n = 0; n < tetEdgeIDSJobArrayExternal.Length; n++)
		{
			tetEdgeIDSJobArrayExternal[n] = tetEdgeIdsVoxelMesh[n];
		}

		for (int n = 0; n < tetSurfaceTriIdsJobArrayExternal.Length; n++)
		{
			tetSurfaceTriIdsJobArrayExternal[n] = tetSurfaceTriIdsVoxelMesh[n];
		}

		//Debug.Log(((Time.realtimeSinceStartup-startTime1)*1000f)+" ms");
		
		combineVoxelsJob2 combineVoxelsJob = new combineVoxelsJob2
		{
			vertsJobArrayInternal = vertsJobArrayExternal,
			tetIDSJobArrayInternal = tetIDSJobArrayExternal,
			tetEdgeIDSJobArrayInternal = tetEdgeIDSJobArrayExternal,
			tetSurfaceTriIdsJobArrayInternal = tetSurfaceTriIdsJobArrayExternal
		};
		//Debug.Log(((Time.realtimeSinceStartup-startTime1)*1000f)+" ms");
		JobHandle jobHandle = combineVoxelsJob.Schedule();
		jobHandle.Complete();
		//Debug.Log(((Time.realtimeSinceStartup-startTime1)*1000f)+" ms");

		for (int n = 0; n < vertsJobArrayExternal.Length; n++)
		{
			vertsVoxelMesh[n] = vertsJobArrayExternal[n];
		}

		for (int n = 0; n < tetIDSJobArrayExternal.Length; n++)
		{
			tetIdsVoxelMesh[n] = tetIDSJobArrayExternal[n];
		}

		for (int n = 0; n < tetEdgeIDSJobArrayExternal.Length; n++)
		{
			tetEdgeIdsVoxelMesh[n] = tetEdgeIDSJobArrayExternal[n];
		}

		for (int n = 0; n < tetSurfaceTriIdsJobArrayExternal.Length; n++)
		{
			tetSurfaceTriIdsVoxelMesh[n] = tetSurfaceTriIdsJobArrayExternal[n];
		}

		vertsJobArrayExternal.Dispose();
		tetIDSJobArrayExternal.Dispose();
		tetEdgeIDSJobArrayExternal.Dispose();
		tetSurfaceTriIdsJobArrayExternal.Dispose();
	}

	[BurstCompile]
	public struct combineVoxelsJob2 : IJob
	{
		public NativeArray<float> vertsJobArrayInternal;
		public NativeArray<int> tetIDSJobArrayInternal;
		public NativeArray<int> tetEdgeIDSJobArrayInternal;
		public NativeArray<int> tetSurfaceTriIdsJobArrayInternal;
		//public NativeHashSet<int> processedIndicesInternal = new NativeHashSet<int>;

		public void Execute()
		{
			for (int i = 0; i < vertsJobArrayInternal.Length/3; i++)
			{
				for (int j = 0; j < vertsJobArrayInternal.Length/3; j++)
				{
					if(vertsJobArrayInternal[3*i]==vertsJobArrayInternal[3*j] && vertsJobArrayInternal[3*i+1]==vertsJobArrayInternal[3*j+1] && vertsJobArrayInternal[3*i+2]==vertsJobArrayInternal[3*j+2] && i!=j)
					{
						//connectionCount++;
						int iPos = 3*i;
						int jPos = 3*j;
						
						for (int k = 0; k < tetEdgeIDSJobArrayInternal.Length; k++)
						{
							if (tetEdgeIDSJobArrayInternal[k] == jPos/3)
							{
								tetEdgeIDSJobArrayInternal[k] = iPos/3;
							}
						}

						for (int l = 0; l < tetIDSJobArrayInternal.Length; l++)
						{
							if (tetIDSJobArrayInternal[l] == jPos/3)
							{
								tetIDSJobArrayInternal[l] = iPos/3;
							}
						}

						for (int m = 0; m < tetSurfaceTriIdsJobArrayInternal.Length; m++)
						{
							if (tetSurfaceTriIdsJobArrayInternal[m] == jPos/3)
							{
								tetSurfaceTriIdsJobArrayInternal[m] = iPos/3;
							}
						}
					}
				}
			}
		}
	}

	[BurstCompile]
	public struct combineVoxelsJob : IJob
	{
		public NativeArray<float> vertsJobArrayInternal;
		public NativeArray<int> tetIDSJobArrayInternal;
		public NativeArray<int> tetEdgeIDSJobArrayInternal;
		public NativeArray<int> tetSurfaceTriIdsJobArrayInternal;
		//public NativeHashSet<int> processedIndicesInternal = new NativeHashSet<int>;

		public void Execute()
		{
			//NativeHashSet<float> processedIndices = new NativeHashSet<float>();
			NativeHashSet<int> processedIndicesInternal = new NativeHashSet<int>();
			Debug.Log("hello");  
			//NativeHashSet<int> processedIndices;
			//HashSet<int> processedIndices = new HashSet<int>();

			for (int i = 0; i < vertsJobArrayInternal.Length / 3; i++)
			{
				if (processedIndicesInternal.Contains(i))
				{
					continue; // Skip already processed indices
				}

				float3 iPosition = new float3(vertsJobArrayInternal[3 * i], vertsJobArrayInternal[3 * i + 1], vertsJobArrayInternal[3 * i + 2]);

				for (int j = i + 1; j < vertsJobArrayInternal.Length / 3; j++)
				{
					float3 jPosition = new float3(vertsJobArrayInternal[3 * j], vertsJobArrayInternal[3 * j + 1], vertsJobArrayInternal[3 * j + 2]);
					if (iPosition.x==jPosition.x && iPosition.y==jPosition.y && iPosition.z==jPosition.z)
					{
						processedIndicesInternal.Add(j);

						int iPos = 3 * i;
						int jPos = 3 * j;
						
						for (int k = 0; k < tetEdgeIDSJobArrayInternal.Length; k++)
						{
							if (tetEdgeIDSJobArrayInternal[k] == jPos / 3)
							{
								tetEdgeIDSJobArrayInternal[k] = iPos / 3;
							}
						}
						
						for (int l = 0; l < tetIDSJobArrayInternal.Length; l++)
						{
							if (tetIDSJobArrayInternal[l] == jPos / 3)
							{
								tetIDSJobArrayInternal[l] = iPos / 3;
							}
						}

						for (int m = 0; m < tetSurfaceTriIdsJobArrayInternal.Length; m++)
						{
							if (tetSurfaceTriIdsJobArrayInternal[m] == jPos / 3)
							{
								tetSurfaceTriIdsJobArrayInternal[m] = iPos / 3;
							}
						}
						
					}
				}
			}
		}
	}

	[BurstCompile]
	public struct tetIDSJob : IJob
	{
		public NativeArray<int> tetIDSJobArrayInternal;
		[ReadOnly] public int i1;
		[ReadOnly] public int j1;

		public void Execute()
		{
			for (int k = 0; k < tetIDSJobArrayInternal.Length; k++)
			{
				if (tetIDSJobArrayInternal[k] == j1 / 3)
				{
					tetIDSJobArrayInternal[k] = i1 / 3;
				}
			}
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

	private void combineVoxels() //Legacy code that has been superseeded by combineAndOptimiseVoxels()
	{
		for (int i = 0; i < vertsVoxelMesh.Length/3; i++)
		{
			for (int j = 0; j < vertsVoxelMesh.Length/3; j++)
			{
				if(vertsVoxelMesh[3*i]==vertsVoxelMesh[3*j] && vertsVoxelMesh[3*i+1]==vertsVoxelMesh[3*j+1] && vertsVoxelMesh[3*i+2]==vertsVoxelMesh[3*j+2] && i!=j)
				{
					connectionCount++;
					int iPos = 3*i;
					int jPos = 3*j;
					
					for (int k = 0; k < tetEdgeIdsVoxelMesh.Length; k++)
					{
						if (tetEdgeIdsVoxelMesh[k] == jPos/3)
						{
							tetEdgeIdsVoxelMesh[k] = iPos/3;
						}
					}

					for (int l = 0; l < tetIdsVoxelMesh.Length; l++)
					{
						if (tetIdsVoxelMesh[l] == jPos/3)
						{
							tetIdsVoxelMesh[l] = iPos/3;
						}
					}

					for (int m = 0; m < tetSurfaceTriIdsVoxelMesh.Length; m++)
					{
						if (tetSurfaceTriIdsVoxelMesh[m] == jPos/3)
						{
							tetSurfaceTriIdsVoxelMesh[m] = iPos/3;
						}
					}
				}
			}
		}
	}
}
