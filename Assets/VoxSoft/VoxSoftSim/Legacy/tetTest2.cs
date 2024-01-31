//Using this test to add a second tetrahedron to the mesh

using System.Collections;
using System.Collections.Generic;
using JetBrains.Annotations;
using UnityEngine;

public class tetTest2 : TetrahedronData
{
	private readonly float[] verts;

	//Getters
	public override float[] GetVerts => verts;

	public override int[] GetTetIds => tetIds;

	public override int[] GetTetEdgeIds => tetEdgeIds;

	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIds;



	public tetTest2()
	{
		//Convert from double to float
		verts = new float[vertsDouble.Length];

		for (int i = 0; i < verts.Length; i++)
		{
			verts[i] = (float)vertsDouble[i];
        }
    }



	//Vertices (x, y, z) come after each other so divide by 3 to get total vertices
	//Provides the vertices of each particle in the mesh
	private double[] vertsDouble =
	{
		1f,0f,0f, 
		-0.5f,Mathf.Sqrt(3)/2f,0f, 
		-0.5f,0,0f,
		0f, 0f, Mathf.Sqrt(2),
		0f,0f,-1f
	};

	//Provides the ID position of the vertices that make up a tetrahedron
	private int[] tetIds =
	{
		0,1,2,3, 0,2,1,4
	};

	//Provides the connections between each one of the edges in a tetrahedron,
	//unlike tetIds the edges should not be repeated with connecting tetrahedrals as they will be looped over.
	private int[] tetEdgeIds =
	{
		//0,1, 1,2, 2,0, 0,3, 2,3, 1,3
		//0,2, 2,1, 1,0, 1,3, 3,0, 2,3,    0,4, 1,4, 2,4, 2,0, 0,1, 1,2 
		0,2, 2,1, 1,0, 1,3, 3,0, 2,3,    0,4, 1,4, 2,4
	};

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		//1,3,2, 0,2,3, 0,3,1, 0,1,2
		0,2,1, 1,2,3, 1,3,0, 0,3,2,    0,2,4, 4,2,1, 4,1,0, 0,1,2
	};
}
