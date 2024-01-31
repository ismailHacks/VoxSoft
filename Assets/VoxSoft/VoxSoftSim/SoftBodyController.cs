using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;


//Simple and unbreakable simulation of soft bodies using Extended Position Based Dynamics (XPBD)
//Based on https://matthias-research.github.io/pages/tenMinutePhysics/index.html
public class SoftBodyController : MonoBehaviour
{
    //Public
    public GameObject softBodyMeshPrefabGO;
    
    public Texture2D cursorTexture;


    //Private
    private readonly List<SoftBodySimulationVectors> allSoftBodies = new ();

    private int numberOfBodies = 1;

    private const int SEED = 0;

    //What we use to grab the balls
    private Grabber grabber;

    private bool simulate = true;

    private readonly Color[] colors = new Color[] { Color.red, Color.blue, Color.green, Color.yellow, Color.cyan };

    

    private void Start()
    {
        Random.InitState(SEED);
        
        TetrahedronData softBodyMesh = new voxelTet();
        //TetrahedronData softBodyMesh = new StanfordBunny();


        for (int i = 0; i < numberOfBodies; i++)
        {
            GameObject meshI = Instantiate(softBodyMeshPrefabGO);

            MeshFilter meshFilter = meshI.GetComponent<MeshFilter>();


            //Random pos
            float halfPlayground = 5f;

            float randomX = Random.Range(-halfPlayground, halfPlayground);
            float randomZ = Random.Range(-halfPlayground, halfPlayground);

            Vector3 startPos = new Vector3(0f, 0f, 0f);


            //Random scale
            //float meshScale = Random.Range(2f, 5f);
            float meshScale = 1f;


            //Random color
            MeshRenderer mr = meshI.GetComponent<MeshRenderer>();

            Material mat = mr.material;

            mat.color = colors[Random.Range(0, colors.Length)];
            

            //SoftBodySimulationTutorial softBodySim = new SoftBodySimulationTutorial(meshFilter, softBodyMesh, startPos, meshScale);
            SoftBodySimulationVectors softBodySim = new SoftBodySimulationVectors(meshFilter, softBodyMesh, startPos, meshScale);


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
            softBody.MyFixedUpdate();
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
