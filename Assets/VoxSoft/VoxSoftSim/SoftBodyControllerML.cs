using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;
using Unity.MLAgents;
using Unity.MLAgents.Actuators;
using Unity.MLAgents.Sensors;
using UnityEngine.UIElements;
using UnityEngine.Analytics;


//Simple and unbreakable simulation of soft bodies using Extended Position Based Dynamics (XPBD)
public class SoftBodyControllerML : Agent
{
    //Public
    public GameObject softBodyMeshPrefabGO;
    GameObject meshI;
    public Texture2D cursorTexture;
    public int numSubSteps;
    private int episodeCount = 0;

	//Soft body behavior settings
	//Compliance (alpha) is the inverse of physical stiffness (k)
	//alpha = 0 means infinitely stiff (hard)
    public float edgeCompliance;
	//Should be 0 or the mesh becomes very flat even for small values 
    public float volCompliance;
    public float dampingCoefficient;

    private float pressure;
    public float pressureControl;
    public float scale;

    //Private
    private readonly List<SoftBodySimulationVectors> allSoftBodies = new ();
    private int numberOfBodies = 1;
    private const int SEED = 0;

    //What we use to grab the particles
    private bool simulate = true;
    private readonly Color[] colors = new Color[] { Color.green, Color.blue, Color.red, Color.yellow, Color.cyan };
    float[] voxPos = new float[343];
    private int stepCount = 0;
    private Vector3 startPos = new Vector3(0,0,0);


    public override void OnEpisodeBegin()
    {
        stepCount = 0;
    }

    public override void OnActionReceived(ActionBuffers actions)
    {
        if(stepCount == 0)
        {
            initialise(actions);
            stepCount++;
        } 
        else
        {
            //Debug.Log("actions out = " + actions.ContinuousActions[voxPos.Length+1]);
            float mappedValue = Mathf.Lerp(-100, 300, Mathf.InverseLerp(-1, 1, actions.ContinuousActions[voxPos.Length+1]));
            pressure = mappedValue;
        }
    }

    public override void CollectObservations(VectorSensor sensor)
    {
        if(stepCount != 0)
        {
            sensor.AddObservation(allSoftBodies[0].CalculateCenterOfMass());
            AddReward(-0.005f);
            float distance = (startPos - allSoftBodies[0].CalculateCenterOfMass()).magnitude;
            AddReward(distance*500);
        }
        //sensor.AddObservation(Position of the thing);
    }

    public override void Heuristic(in ActionBuffers actionsOut)
    {
        ActionSegment<float> continuousActions = actionsOut.ContinuousActions;
        continuousActions[voxPos.Length+1] = pressureControl;
    }

    private void Update()
    {
        if(stepCount != 0)
        {
            foreach (SoftBodySimulationVectors softBody in allSoftBodies)
            {
                softBody.MyUpdate();
            }
        }
    }

    private void FixedUpdate()
    {
        if(stepCount != 0)
        {
            Time.timeScale=0.45f;
            if (!simulate)
            {
                return;
            }

            foreach (SoftBodySimulationVectors softBody in allSoftBodies)
            {
                softBody.MyFixedUpdate(numSubSteps,edgeCompliance, volCompliance, dampingCoefficient, pressure);
                Debug.DrawRay(softBody.CalculateCenterOfMass(), new Vector3(0,1,0) , Color.blue);
            }
        }
    }

    private void initialise(ActionBuffers actions)
    {
        // Destroy the previous mesh if it exists
        if (meshI != null)
        {
            Destroy(meshI);
        }

        // Destroy any existing soft body simulations and clear the list
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            Mesh mesh = softBody.MyOnDestroy();
            Destroy(mesh);
        }
        allSoftBodies.Clear();

        for (int i = 0; i < voxPos.Length; i++)
        {
            voxPos[i] = actions.ContinuousActions[i];
        }

        // Reinitialize the random seed
        UnityEngine.Random.InitState(SEED);
        TetrahedronData softBodyMesh = new voxelTet(scale, voxPos);
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
        startPos = allSoftBodies[0].CalculateCenterOfMass();
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
