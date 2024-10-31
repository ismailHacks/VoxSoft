using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;


//Simple and unbreakable simulation of soft bodies using Extended Position Based Dynamics (XPBD)
public class SoftBodyControllerA : MonoBehaviour
{
    //Public
    public GameObject softBodyMeshPrefabGO;
    GameObject meshI;
    public Texture2D cursorTexture;
    public int numSubSteps;

	//Soft body behavior settings
	//Compliance (alpha) is the inverse of physical stiffness (k)
	//alpha = 0 means infinitely stiff (hard)
    public float edgeCompliance;
    public float volCompliance;
    public float dampingCoefficient;

    public float pressure;
    public float scale;

    //Private
    private readonly List<SoftBodySimulationVectors> allSoftBodies = new ();
    private int numberOfBodies = 1;
    private const int SEED = 0;

    //What we use to grab the particles
    private Grabber grabber;
    private bool simulate = true;
    private readonly Color[] colors = new Color[] { Color.green, Color.blue, Color.red, Color.yellow, Color.cyan };
    float[] voxPos = new float[343];
    private int counter = 0;
    private int counterMax = 1000;


    private string exportWidth = "/Data/Width";
    private string exportHeight = "/Data/Heigh";
    private string exportWallThickness = "/Data/WallThickness";
    private string exportFitness = "/Data/Fitness";

    private void Start()
    {
        UnityEngine.Random.InitState(SEED);
        
        TetrahedronData softBodyMesh = new voxelTet(scale, voxPos, 2, 2, 2);
        voxelTet voxelSpecific = (voxelTet)softBodyMesh;
        //TetrahedronData softBodyMesh = new tetGen(scale);
        //TetrahedronData softBodyMesh = new StanfordBunny();

        for (int i = 0; i < numberOfBodies; i++)
        {
            meshI = Instantiate(softBodyMeshPrefabGO);

            MeshFilter meshFilter = meshI.GetComponent<MeshFilter>();

            Vector3 startPos = new Vector3(0f, 0f, 0f);

            //Random color
            MeshRenderer mr = meshI.GetComponent<MeshRenderer>();

            Material mat = mr.material;

            mat.color = colors[0];
            
            SoftBodySimulationVectors softBodySim = new SoftBodySimulationVectors(meshFilter, softBodyMesh, startPos, scale, voxelSpecific);

            allSoftBodies.Add(softBodySim);
        }

        //Init the grabber
        grabber = new Grabber(Camera.main);

        Cursor.visible = true;

        Cursor.SetCursor(cursorTexture, Vector2.zero, CursorMode.ForceSoftware);
    }

    private void reinitialise(int width, int height, int wallThickness)
    {
        Debug.Log("reinitialising");

        if (meshI != null)
        {
            Debug.Log("Destroying Mesh");
            Destroy(meshI);
        }
        // Destroy any existing soft body simulations and clear the list
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            Mesh mesh = softBody.MyOnDestroy();
            Destroy(mesh);
        }

        allSoftBodies.Clear();

        //Reinitialize the random seed
        //UnityEngine.Random.InitState(SEED);
        TetrahedronData softBodyMesh = new voxelTet(scale, voxPos, width, height, wallThickness);
        voxelTet voxelSpecific = (voxelTet)softBodyMesh;

        // Instantiate a new soft body mesh
        for (int i = 0; i < numberOfBodies; i++)
        {
            meshI = Instantiate(softBodyMeshPrefabGO);

            MeshFilter meshFilter = meshI.GetComponent<MeshFilter>();

            Vector3 startPos = new Vector3(0f, 0f, 0f);

            // Set random color
            MeshRenderer mr = meshI.GetComponent<MeshRenderer>();

            Material mat = mr.material;

            mat.color = colors[0];

            SoftBodySimulationVectors softBodySim = new SoftBodySimulationVectors(meshFilter, softBodyMesh, startPos, scale, voxelSpecific);

            allSoftBodies.Add(softBodySim);
        }


    }

    private void Update()
    {
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            softBody.MyUpdate();
        }

        grabber.MoveGrab();

        //Pause simulation
        if (Input.GetKeyDown(KeyCode.P))
        {
            simulate = !simulate;
        }
    }

    private void LateUpdate()
    {
        if (Input.GetMouseButtonDown(0))
        {
            List<IGrabbable> temp = new List<IGrabbable>(allSoftBodies);
        
            grabber.StartGrab(temp);
        }

        if (Input.GetMouseButtonUp(0))
        {
            grabber.EndGrab();
        }
    }


    private int i = 0;
    private int j = 0;
    private int k = 0;

    private void FixedUpdate()
    {
        if (counter > counterMax)
        {
            addRecord(i+2, (Application.dataPath + exportWidth));
            addRecord(j+2, (Application.dataPath + exportHeight));
            addRecord(k+1, (Application.dataPath + exportWallThickness));
            addRecord(allSoftBodies[0].GetAverageTopVoxelsVerticalDisplacement(), (Application.dataPath + exportFitness));

            reinitialise(i + 2, j + 2, k + 1); //width, height, Wall Thickness
            counter = 0;

            // Increment k and handle carry-over
            k++;
            if (k >= 3)
            {
                k = 0;
                j++;
                if (j >= 6)
                {
                    j = 0;
                    i++;
                    if (i >= 9)
                    {
                        // Reset or handle completion
                        //i = 0;
                        // Optionally, set simulate = false if you want to stop
                        simulate = false;
                    }
                }
            }
        }

        if (!simulate)
        {
            return;
        }

        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            softBody.MyFixedUpdate(numSubSteps, edgeCompliance, volCompliance, dampingCoefficient, pressure);
        }
        counter++;
    }

    private void OnDestroy()
    {
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            Mesh mesh = softBody.MyOnDestroy();
            Destroy(mesh);
        }
    }

    public static void addRecord(float fitness, string filepath)
    {
        try
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath, true))
            {
                string str = fitness.ToString();
                file.WriteLine(str);
            }
        }
        catch(Exception ex)
        {
            throw new ApplicationException("error :", ex);
        }
    }
}
