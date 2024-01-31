using UnityEngine;

[RequireComponent(typeof(MeshFilter))]
public class TetrahedronMeshGenerator : MonoBehaviour
{
    void Start()
    {
        GenerateTetrahedronMesh();
    }

    void GenerateTetrahedronMesh()
    {
        MeshFilter meshFilter = GetComponent<MeshFilter>();
        Mesh mesh = meshFilter.mesh;
        mesh.Clear();

        Vector3[] vertices = new Vector3[4];

        // Define the vertices of the tetrahedron
        vertices[0] = new Vector3(1f, 0f, 0f);
        vertices[1] = new Vector3(-0.5f, Mathf.Sqrt(3) / 2f, 0f);
        vertices[2] = new Vector3(-0.5f, -Mathf.Sqrt(3) / 2f, 0f);
        vertices[3] = new Vector3(0f, 0f, Mathf.Sqrt(2));
        /*int[] triangles = new int[12];

        // Define the triangles that make up the tetrahedron
        triangles[0] = 0;
        triangles[1] = 1;
        triangles[2] = 2;

        triangles[3] = 0;
        triangles[4] = 2;
        triangles[5] = 3;

        triangles[6] = 0;
        triangles[7] = 3;
        triangles[8] = 1;

        triangles[9] = 1;
        triangles[10] = 3;
        triangles[11] = 2;*/

        // Assign vertices and triangles to the mesh
        double[] vertsDouble =
        {
            1f, 0f, 0f, -0.5f, Mathf.Sqrt(3) / 2f, 0f, -0.5f, -Mathf.Sqrt(3) / 2f, 0f, 0f, 0f, Mathf.Sqrt(2)
        };

        int[] tetSurfaceTriIds =
        {
            0,1,2, 0,2,3, 0,3,1, 1,3,2
        };

        mesh.vertices = vertices;
        mesh.triangles = tetSurfaceTriIds;

        // Automatically calculate normals and bounds
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
    }
}

