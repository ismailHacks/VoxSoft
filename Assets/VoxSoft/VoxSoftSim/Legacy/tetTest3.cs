//Using this test to add tetrahedrons in a cylinder like pattern

using System.Collections;
using System.Collections.Generic;
using JetBrains.Annotations;
using Unity.Mathematics;
using UnityEngine;

public class tetTest3 : TetrahedronData
{
	private readonly float[] verts;

	//Getters
	public override float[] GetVerts => vertsCylinder;
	public override int[] GetTetIds => tetIdsCylinder;
	public override int[] GetTetEdgeIds => tetEdgeIdsCylinder;
	public override int[] GetTetSurfaceTriIds => tetSurfaceTriIdsCylinder;

	public float cylinderOuterRadius = 1f;
	public float cylinderInnerRadius = 0.5f;

	private static int rotationNo = 7;
	private float[] vertsCylinder = new float[12*rotationNo];
	private int[] tetIdsCylinder = new int[4*rotationNo];
	private int[] tetEdgeIdsCylinder = new int[12*rotationNo];
	private int[] tetSurfaceTriIdsCylinder = new int[12*rotationNo];




	public tetTest3()
	{
		float angleBetweenTets = (2*(cylinderOuterRadius-cylinderInnerRadius)*Mathf.Tan(Mathf.Deg2Rad*30))/cylinderOuterRadius;
		angleBetweenTets = 360/(Mathf.Rad2Deg*angleBetweenTets);
		float angleBetweenTetsMod = 360/(int)angleBetweenTets;
		angleBetweenTetsMod = Mathf.Deg2Rad*angleBetweenTetsMod;
		
		for (int i = 0; i < rotationNo; i++)
		{
			float currentAngle = angleBetweenTetsMod*i;
			if (i==0)
			{
				vertsCylinder[0] = cylinderOuterRadius*Mathf.Sin(currentAngle);
				vertsCylinder[1] = 0;
				vertsCylinder[2] = cylinderOuterRadius*Mathf.Cos(currentAngle);
				vertsCylinder[3] = cylinderInnerRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[4] = 0;
				vertsCylinder[5] = cylinderInnerRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[6] = cylinderOuterRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[7] =0.5f;
				vertsCylinder[8] = cylinderOuterRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[9] = cylinderOuterRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod);
				vertsCylinder[10] = 0;
				vertsCylinder[11] = cylinderOuterRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod);

				tetIdsCylinder[0] = i;
				tetIdsCylinder[1] = i+1;
				tetIdsCylinder[2] = i+2;
				tetIdsCylinder[3] = i+3;

				tetEdgeIdsCylinder[0] = i;
				tetEdgeIdsCylinder[1] = i+1;
				tetEdgeIdsCylinder[2] = i;
				tetEdgeIdsCylinder[3] = i+2;
				tetEdgeIdsCylinder[4] = i;
				tetEdgeIdsCylinder[5] = i+3;
				tetEdgeIdsCylinder[6] = i+1;
				tetEdgeIdsCylinder[7] = i+2;
				tetEdgeIdsCylinder[8] = i+1;
				tetEdgeIdsCylinder[9] = i+3;
				tetEdgeIdsCylinder[10] = i+2;
				tetEdgeIdsCylinder[11] = i+3;

				tetSurfaceTriIdsCylinder[0] = i;
				tetSurfaceTriIdsCylinder[1] = i+2;
				tetSurfaceTriIdsCylinder[2] = i+1;
				tetSurfaceTriIdsCylinder[3] = i+1;
				tetSurfaceTriIdsCylinder[4] = i+2;
				tetSurfaceTriIdsCylinder[5] = i+3;
				tetSurfaceTriIdsCylinder[6] = i+1;
				tetSurfaceTriIdsCylinder[7] = i+3;
				tetSurfaceTriIdsCylinder[8] = i+0;
				tetSurfaceTriIdsCylinder[9] = i+0;
				tetSurfaceTriIdsCylinder[10] = i+3;
				tetSurfaceTriIdsCylinder[11] = i+2;
			}
			else
			{
				vertsCylinder[9*i+3] = cylinderInnerRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[9*i+4] = 0;
				vertsCylinder[9*i+5] = cylinderInnerRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[9*i+6] = cylinderOuterRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[9*i+7] =0.5f;
				vertsCylinder[9*i+8] = cylinderOuterRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod/2);
				vertsCylinder[9*i+9] = cylinderOuterRadius*Mathf.Sin(currentAngle+angleBetweenTetsMod);
				vertsCylinder[9*i+10] = 0;
				vertsCylinder[9*i+11] = cylinderOuterRadius*Mathf.Cos(currentAngle+angleBetweenTetsMod);
				
				tetIdsCylinder[4*i+0] = 3*i;
				tetIdsCylinder[4*i+1] = 3*i+1;
				tetIdsCylinder[4*i+2] = 3*i+2;
				tetIdsCylinder[4*i+3] = 3*i+3;

				tetEdgeIdsCylinder[12*i+0] = 3*i;
				tetEdgeIdsCylinder[12*i+1] = 3*i+1;
				tetEdgeIdsCylinder[12*i+2] = 3*i;
				tetEdgeIdsCylinder[12*i+3] = 3*i+2;
				tetEdgeIdsCylinder[12*i+4] = 3*i;
				tetEdgeIdsCylinder[12*i+5] = 3*i+3;
				tetEdgeIdsCylinder[12*i+6] = 3*i+1;
				tetEdgeIdsCylinder[12*i+7] = 3*i+2;
				tetEdgeIdsCylinder[12*i+8] = 3*i+1;
				tetEdgeIdsCylinder[12*i+9] = 3*i+3;
				tetEdgeIdsCylinder[12*i+10] = 3*i+2;
				tetEdgeIdsCylinder[12*i+11] = 3*i+3;

				tetSurfaceTriIdsCylinder[12*i+0] = 3*i+0;
				tetSurfaceTriIdsCylinder[12*i+1] = 3*i+2;
				tetSurfaceTriIdsCylinder[12*i+2] = 3*i+1;
				tetSurfaceTriIdsCylinder[12*i+3] = 3*i+1;
				tetSurfaceTriIdsCylinder[12*i+4] = 3*i+2;
				tetSurfaceTriIdsCylinder[12*i+5] = 3*i+3;
				tetSurfaceTriIdsCylinder[12*i+6] = 3*i+1;
				tetSurfaceTriIdsCylinder[12*i+7] = 3*i+3;
				tetSurfaceTriIdsCylinder[12*i+8] = 3*i+0;
				tetSurfaceTriIdsCylinder[12*i+9] = 3*i+0;
				tetSurfaceTriIdsCylinder[12*i+10] = 3*i+3;
				tetSurfaceTriIdsCylinder[12*i+11] = 3*i+2;
				Debug.Log(vertsCylinder);
			}
		}
	}

	//Provides the ID position of the vertices that make up a tetrahedron
	/*private int[] tetIds =
	{
		0,1,2,3, 3,4,5,6
	};*/

	//Provides the connections between each one of the edges in a tetrahedron,
	//unlike tetIds the edges should not be repeated with connecting tetrahedrals as they will be looped over.
	/*private int[] tetEdgeIds =
	{
		0,1, 0,2, 0,3, 1,2, 1,3, 2,3,   3,4, 3,5, 3,6, 4,5, 4,6, 5,6
	};*/

	//Provides the connections between all surfaces that are visible in order to render a mesh, must be clockwise done when looking at the surface.
	private int[] tetSurfaceTriIds =
	{
		0,2,1, 1,2,3, 1,3,0, 0,3,2,    3,5,4, 4,5,6, 4,6,3, 3,6,5
	};
}
