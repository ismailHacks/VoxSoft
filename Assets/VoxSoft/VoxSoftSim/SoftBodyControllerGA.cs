using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;
using UserInteraction;


//Simple and unbreakable simulation of soft bodies using Extended Position Based Dynamics (XPBD)
public class SoftBodyControllerGA : MonoBehaviour
{
    //Public
    public GameObject softBodyMeshPrefabGO;
    public Texture2D cursorTexture;
    public int numSubSteps;

    [Header("Genetic Algorithm")]
    static int populationSize = 30;
    private float mutationRate = 0.05f;
    private int elitism = 2;
    private GeneticAlgorithm<float> ga;
    private System.Random random;
    private int DNAno = 0;
    private float[] outputGenesEdgeC = new float[populationSize];
    private float[] outputGenesVolC = new float[populationSize];
    private float[] outputGenesStep = new float[populationSize];
    private float[] Fitness = new float[populationSize];
    private string exportEdgeC = "/Data/EdgeCompliance";
    private string exportVolC = "/Data/VolumeCompliance";
    private string exportStep = "/Data/Steps";
    private string exportFitness = "/Data/Fitness";
    private string exportMatchFitness = "/Data/MatchFitness";

    public float[] matchFitness = new float[populationSize];

    int counter = 0;
    private int counterLimit = 3000;

    private float currentBestFitness = 0f;
    private int generation = 0;

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
            
            SoftBodySimulationVectors softBodySim = new SoftBodySimulationVectors(meshFilter, softBodyMesh, startPos);

            allSoftBodies.Add(softBodySim);
        }

        //Init the grabber
        grabber = new Grabber(Camera.main);

        Cursor.visible = true;

        Cursor.SetCursor(cursorTexture, Vector2.zero, CursorMode.ForceSoftware);

        random = new System.Random();
        ga = new GeneticAlgorithm<float>(populationSize, 2, random, getRandomNo, FitnessFunction, elitism, mutationRate);
        //ga = new GeneticAlgorithm<float>(populationSize, 3, random, getRandomNo, FitnessFunction, elitism, mutationRate);


        for (int i = 0; i < populationSize; i++) 
        {
            outputGenesEdgeC[i] = scaleVal(0f,1f,0f,0.3f, ga.Population[i].Genes[0]);
            outputGenesVolC[i] = scaleVal(0f,1f,0f,0.1f, ga.Population[i].Genes[1]);
            //outputGenesStep[i] = scaleVal(0f,1f,150f,150f, ga.Population[i].Genes[2]);
        }
        edgeCompliance = outputGenesEdgeC[0];
        volCompliance = outputGenesVolC[0];
        //numSubSteps = (int)outputGenesStep[0];
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
        counter++;
        if(allSoftBodies[0].converged || counter > counterLimit)
        {
            matchFitness[DNAno] = allSoftBodies[0].fitnessCalculate();
            counter = 0; 
            DNAno++;
            if(DNAno > populationSize-1)
            {
                DNAno = 0;
                addRecord(ga.BestFitness, (Application.dataPath + exportFitness));
                addRecordGA(matchFitness, (Application.dataPath  + exportMatchFitness));
                addRecordGA(outputGenesEdgeC, (Application.dataPath  + exportEdgeC));
                addRecordGA(outputGenesVolC, (Application.dataPath  + exportVolC));
                //addRecordGA(outputGenesStep, (Application.dataPath  + exportStep));

                ga.NewGeneration();
                
                for (int i = 0; i < populationSize; i++) 
                {
                    outputGenesEdgeC[i] = scaleVal(0f,1f,0f,0.3f, ga.Population[i].Genes[0]);
                    outputGenesVolC[i] = scaleVal(0f,1f,0f,0.1f, ga.Population[i].Genes[1]);
                    //outputGenesStep[i] = scaleVal(0f,1f,10f,100f, ga.Population[i].Genes[2]);
                }
                
                edgeCompliance = outputGenesEdgeC[DNAno];
                volCompliance = outputGenesVolC[DNAno];
                //numSubSteps = (int)outputGenesStep[DNAno];
            }
             else
            {
                edgeCompliance = outputGenesEdgeC[DNAno];
                volCompliance = outputGenesVolC[DNAno];
                //numSubSteps = (int)outputGenesStep[DNAno];
            }
        }
    }

    private void OnDestroy()
    {
        foreach (SoftBodySimulationVectors softBody in allSoftBodies)
        {
            Mesh mesh = softBody.MyOnDestroy();

            Destroy(mesh);
        }
    }
    
    private float FitnessFunction(int index)
    {
        float score = 0;
        //score = (1-Mathf.Abs((displacementFitness[index]-targetDisplacement)))+(1-stressFitness[index])/20f; //+1-stress induced
        score = matchFitness[index]; //+1-stress induced
        return score;
    }

    private float getRandomNo()
    {
        float number = UnityEngine.Random.Range(0f, 1f);
        return number;
    } 

    public float scaleVal(float OldMin, float OldMax, float NewMin, float NewMax, float OldValue)
    {
        float OldRange = (OldMax - OldMin);
        float NewRange = (NewMax - NewMin);
        float NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin;
        return(NewValue);
    }

    public static float Sigmoid(float value)
    {
        return 1.0f / (1.0f + Mathf.Exp(-value));
    }

    public static void addRecordGA(float[] variable ,string filepath)
    {
        try
        {
            using (System.IO.StreamWriter file = new System.IO.StreamWriter(@filepath, true))
            {
                string str = null;
                for (int i = 0; i < populationSize; i++) 
                {
                    string concat = String.Concat(str, variable[i].ToString(), ",");
                    str = concat;
                }
                file.WriteLine(str);
            }
        }
        catch(Exception ex)
        {
            throw new ApplicationException("error :", ex);
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
