using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;
using Unity.MLAgents;
using Unity.MLAgents.Actuators;
using Unity.MLAgents.Sensors;
using UnityEngine.UIElements;


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

    public float pressure;
    public float scale;

    //Private
    private readonly List<SoftBodySimulationVectors> allSoftBodies = new ();
    private int numberOfBodies = 1;
    private const int SEED = 0;

    //What we use to grab the particles
    private bool simulate = true;
    private readonly Color[] colors = new Color[] { Color.green, Color.blue, Color.red, Color.yellow, Color.cyan };


    public override void OnEpisodeBegin()
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

        // Reinitialize the random seed
        UnityEngine.Random.InitState(SEED);
        TetrahedronData softBodyMesh = new voxelTet(scale);
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

    public override void OnActionReceived(ActionBuffers actions)
    {
        //Debug.Log(actions.ContinuousActions[0]);
        pressure = actions.ContinuousActions[0];
        //Need to add in displacement of the soft robot
    }

    public override void CollectObservations(VectorSensor sensor)
    {
        //sensor.AddObservation(Position of the thing);
        //SetReward(1f);
    }

    private void Update()
    {
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            softBody.MyUpdate();
        }
    }

    private void FixedUpdate()
    {
        Time.timeScale=0.15f;
        //Timers.Reset();
        if (!simulate)
        {
            return;
        }

        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            softBody.MyFixedUpdate(numSubSteps,edgeCompliance, volCompliance, dampingCoefficient, pressure);
            Debug.DrawRay(softBody.CalculateCenterOfMass(), new Vector3(0,1,0) , Color.blue);
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
