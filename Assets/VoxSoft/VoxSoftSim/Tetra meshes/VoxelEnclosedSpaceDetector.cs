using System.Collections.Generic;
using UnityEngine;

public class VoxelEnclosedSpaceDetector : MonoBehaviour
{
    // Dimensions of the voxel grid
    int sizeX, sizeY, sizeZ;

    // The voxel grid
    bool[,,] voxelData;

    // Visited markers for flood fill
    bool[,,] visited;

    bool debug = false;

    Dictionary<Vector3Int, int> voxelPositionToID;

    public Dictionary<string, List<int>> DetectEnclosedSpaces(
        bool[,,] inputVoxelData, Dictionary<Vector3Int, int> positionToIDMap)
    {
        voxelData = inputVoxelData;
        voxelPositionToID = positionToIDMap;
        sizeX = voxelData.GetLength(0);
        sizeY = voxelData.GetLength(1);
        sizeZ = voxelData.GetLength(2);

        visited = new bool[sizeX, sizeY, sizeZ];

        // Step 1: Flood fill from the outside
        FloodFillFromOutside();

        // Step 2: Find enclosed empty voxels
        List<Vector3Int> enclosedEmptyVoxels = FindEnclosedEmptyVoxels();

        // Step 3: Identify wall voxels and exposed faces
        Dictionary<Vector3Int, List<string>> wallVoxels = FindWallVoxels(enclosedEmptyVoxels);

        // Step 4: Map wall voxels to IDs and group by face direction
        Dictionary<string, List<int>> faceDirectionToVoxelIDs = MapWallVoxelsToIDs(wallVoxels);

        return faceDirectionToVoxelIDs;
    }

    Dictionary<string, List<int>> MapWallVoxelsToIDs(Dictionary<Vector3Int, List<string>> wallVoxels)
    {
        Dictionary<string, List<int>> faceDirectionToVoxelIDs = new Dictionary<string, List<int>>();

        foreach (var kvp in wallVoxels)
        {
            Vector3Int pos = kvp.Key;
            List<string> faces = kvp.Value;

            // Get the voxel ID from the position
            if (voxelPositionToID.TryGetValue(pos, out int voxelID))
            {
                foreach (string face in faces)
                {
                    if (!faceDirectionToVoxelIDs.ContainsKey(face))
                    {
                        faceDirectionToVoxelIDs[face] = new List<int>();
                    }
                    faceDirectionToVoxelIDs[face].Add(voxelID);
                }
            }
            else
            {
                Debug.LogWarning($"Voxel at position {pos} does not have a corresponding ID.");
            }
        }

        return faceDirectionToVoxelIDs;
    }

    // Directions for flood fill (including diagonals)
    Vector3Int[] floodFillDirections = new Vector3Int[]
    {
        new Vector3Int(1, 0, 0),   // Right
        new Vector3Int(-1, 0, 0),  // Left
        new Vector3Int(0, 1, 0),   // Up
        new Vector3Int(0, -1, 0),  // Down
        new Vector3Int(0, 0, 1),   // Forward
        new Vector3Int(0, 0, -1),  // Backward
        new Vector3Int(1, 1, 1),
        new Vector3Int(1, 1, -1),
        new Vector3Int(-1, 1, 1),
        new Vector3Int(-1, 1, -1),
        new Vector3Int(1, -1, 1),
        new Vector3Int(1, -1, -1),
        new Vector3Int(-1, -1, 1),
        new Vector3Int(-1, -1, -1)
    };

    // Directions for neighbor checking (only face directions)
    Vector3Int[] faceDirections = new Vector3Int[]
    {
        new Vector3Int(1, 0, 0),   // Right
        new Vector3Int(-1, 0, 0),  // Left
        new Vector3Int(0, 1, 0),   // Up
        new Vector3Int(0, -1, 0),  // Down
        new Vector3Int(0, 0, 1),   // Forward
        new Vector3Int(0, 0, -1)   // Backward
    };

    // Exposed face names corresponding to face directions
    string[] faceNames = new string[]
    {
        "Right",
        "Left",
        "Top",
        "Bottom",
        "Front",
        "Back"
    };

    void FloodFillFromOutside()
    {
        Queue<Vector3Int> queue = new Queue<Vector3Int>();

        // Enqueue all border empty voxels
        for (int x = 0; x < sizeX; x++)
        {
            for (int y = 0; y < sizeY; y++)
            {
                EnqueueIfEmpty(queue, x, y, 0);
                EnqueueIfEmpty(queue, x, y, sizeZ - 1);
            }
        }
        for (int x = 0; x < sizeX; x++)
        {
            for (int z = 0; z < sizeZ; z++)
            {
                EnqueueIfEmpty(queue, x, 0, z);
                EnqueueIfEmpty(queue, x, sizeY - 1, z);
            }
        }
        for (int y = 0; y < sizeY; y++)
        {
            for (int z = 0; z < sizeZ; z++)
            {
                EnqueueIfEmpty(queue, 0, y, z);
                EnqueueIfEmpty(queue, sizeX - 1, y, z);
            }
        }

        // Flood fill algorithm
        while (queue.Count > 0)
        {
            Vector3Int pos = queue.Dequeue();

            foreach (Vector3Int dir in floodFillDirections)
            {
                Vector3Int neighbor = pos + dir;

                if (IsWithinBounds(neighbor) && !visited[neighbor.x, neighbor.y, neighbor.z] && !voxelData[neighbor.x, neighbor.y, neighbor.z])
                {
                    visited[neighbor.x, neighbor.y, neighbor.z] = true;
                    if(debug)
                    {
                            GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                            cube.transform.position = new Vector3(neighbor.x*0.005f+0.0025f, neighbor.y*0.005f+0.0025f, neighbor.z*0.005f+0.0025f);
                            cube.transform.localScale = new Vector3(0.0025f,0.0025f,0.0025f);
                            Renderer cubeRenderer = cube.GetComponent<Renderer>();
                            cubeRenderer.material.color = Color.red;
                    }
                    queue.Enqueue(neighbor);
                }
            }
        }
    }

    void EnqueueIfEmpty(Queue<Vector3Int> queue, int x, int y, int z)
    {
        if (!voxelData[x, y, z] && !visited[x, y, z])
        {
            visited[x, y, z] = true;
            queue.Enqueue(new Vector3Int(x, y, z));
            if(debug)
            {
                GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                cube.transform.position = new Vector3(x*0.005f+0.0025f, y*0.005f+0.0025f, z*0.005f+0.0025f);
                cube.transform.localScale = new Vector3(0.0025f,0.0025f,0.0025f);
                Renderer cubeRenderer = cube.GetComponent<Renderer>();
                cubeRenderer.material.color = Color.red;
            }
        }
    }

    List<Vector3Int> FindEnclosedEmptyVoxels()
    {
        List<Vector3Int> enclosedEmptyVoxels = new List<Vector3Int>();

        for (int x = 0; x < sizeX; x++)
        {
            for (int y = 0; y < sizeY; y++)
            {
                for (int z = 0; z < sizeZ; z++)
                {
                    if (!voxelData[x, y, z] && !visited[x, y, z])
                    {
                        enclosedEmptyVoxels.Add(new Vector3Int(x, y, z));
                        if(debug)
                        {
                            GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
                            cube.transform.position = new Vector3(x*0.005f+0.0025f, y*0.005f+0.0025f, z*0.005f+0.0025f);
                            cube.transform.localScale = new Vector3(0.0025f,0.0025f,0.0025f);
                            Renderer cubeRenderer = cube.GetComponent<Renderer>();
                            cubeRenderer.material.color = Color.blue;
                        }
                    }
                }
            }
        }

        return enclosedEmptyVoxels;
    }

    Dictionary<Vector3Int, List<string>> FindWallVoxels(List<Vector3Int> enclosedEmptyVoxels)
    {
        Dictionary<Vector3Int, List<string>> wallVoxels = new Dictionary<Vector3Int, List<string>>();

        foreach (Vector3Int emptyPos in enclosedEmptyVoxels)
        {
            for (int i = 0; i < faceDirections.Length; i++)
            {
                Vector3Int dir = faceDirections[i];
                Vector3Int neighborPos = emptyPos + dir;

                if (IsWithinBounds(neighborPos) && voxelData[neighborPos.x, neighborPos.y, neighborPos.z])
                {
                    if (!wallVoxels.ContainsKey(neighborPos))
                    {
                        wallVoxels[neighborPos] = new List<string>();
                    }
                    wallVoxels[neighborPos].Add(faceNames[i]);
                }
            }
        }

        return wallVoxels;
    }

    void OutputWallVoxels(Dictionary<Vector3Int, List<string>> wallVoxels)
    {
        foreach (var kvp in wallVoxels)
        {
            Vector3Int pos = kvp.Key;
            List<string> faces = kvp.Value;

            Debug.Log($"Voxel at position ({pos.x}, {pos.y}, {pos.z}) has exposed faces: {string.Join(", ", faces)}");
        }
    }

    bool IsWithinBounds(Vector3Int pos)
    {
        return pos.x >= 0 && pos.x < sizeX && pos.y >= 0 && pos.y < sizeY && pos.z >= 0 && pos.z < sizeZ;
    }
}