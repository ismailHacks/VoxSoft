using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;


//Simple and unbreakable simulation of soft bodies using Extended Position Based Dynamics (XPBD)
public class SoftBodyController : MonoBehaviour
{
    //Public
    public GameObject softBodyMeshPrefabGO;
    public Texture2D cursorTexture;
    public int numSubSteps;

	//Soft body behavior settings
	//Compliance (alpha) is the inverse of physical stiffness (k)
	//alpha = 0 means infinitely stiff (hard)
    public float edgeCompliance;
	//Should be 0 or the mesh becomes very flat even for small values 
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

    private void Start()
    {
        UnityEngine.Random.InitState(SEED);
        
        TetrahedronData softBodyMesh = new voxelTet(scale);
        //TetrahedronData softBodyMesh = new tetGen(scale);
        //TetrahedronData softBodyMesh = new StanfordBunny();

        for (int i = 0; i < numberOfBodies; i++)
        {
            GameObject meshI = Instantiate(softBodyMeshPrefabGO);

            MeshFilter meshFilter = meshI.GetComponent<MeshFilter>();

            Vector3 startPos = new Vector3(0f, 0f, 0f);

            //Random color
            MeshRenderer mr = meshI.GetComponent<MeshRenderer>();

            Material mat = mr.material;

            mat.color = colors[0];
            
            SoftBodySimulationVectors softBodySim = new SoftBodySimulationVectors(meshFilter, softBodyMesh, startPos, scale);

            allSoftBodies.Add(softBodySim);
        }

        //Init the grabber
        grabber = new Grabber(Camera.main);

        Cursor.visible = true;

        Cursor.SetCursor(cursorTexture, Vector2.zero, CursorMode.ForceSoftware);
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

    private void FixedUpdate()
    {
        //Timers.Reset();
        if (!simulate)
        {
            return;
        }

        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            softBody.MyFixedUpdate(numSubSteps,edgeCompliance, volCompliance, dampingCoefficient, pressure);
        }
        //Timers.Display();
       
    }

    private void OnDestroy()
    {
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            Mesh mesh = softBody.MyOnDestroy();
            Destroy(mesh);
        }
    }
}
