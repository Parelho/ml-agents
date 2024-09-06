using System;
using UnityEngine;
using Unity.MLAgents;
using Unity.MLAgents.Actuators;
using Unity.MLAgentsExamples;
using Unity.MLAgents.Sensors;
using BodyPart = Unity.MLAgentsExamples.BodyPart;
using Random = UnityEngine.Random;

public class WalkerAgent : Agent
{
    [Header("Walk Speed")]
    [Range(1, 10)]
    [SerializeField]
    //The walking speed to try and achieve
    private float m_TargetWalkingSpeed = 10;

    public float MTargetWalkingSpeed // property
    {
        get { return m_TargetWalkingSpeed; }
        set { m_TargetWalkingSpeed = Mathf.Clamp(value, 1f, m_maxWalkingSpeed); }
    }

    const float m_maxWalkingSpeed = 10; //The max walking speed

    //Should the agent sample a new goal velocity each episode?
    //If true, walkSpeed will be randomly set between zero and m_maxWalkingSpeed in OnEpisodeBegin() 
    //If false, the goal velocity will be walkingSpeed
    public bool randomizeWalkSpeedEachEpisode;
    [Header("Body Parts")] public Transform hips;
    public Transform thighL;
    public Transform shinL;
    public Transform footL;
    public Transform thighR;
    public Transform shinR;
    public Transform footR;
    public Transform targetTransform;

    private bool lastFootWasLeft;
    private bool firstFootTouch;

    private float rayRange = 23f;

    //This will be used as a stabilized model space reference point for observations
    //Because ragdolls can move erratically during training, using a stabilized reference transform improves learning

    //The indicator graphic gameobject that points towards the target
    JointDriveController m_JdController;
    EnvironmentParameters m_ResetParams;

    public override void Initialize()
    {
        //Setup each body part
        m_JdController = GetComponent<JointDriveController>();
        m_JdController.SetupBodyPart(hips);
        m_JdController.SetupBodyPart(thighL);
        m_JdController.SetupBodyPart(shinL);
        m_JdController.SetupBodyPart(footL);
        m_JdController.SetupBodyPart(thighR);
        m_JdController.SetupBodyPart(shinR);
        m_JdController.SetupBodyPart(footR);

        targetTransform = transform.Find("Target");

        m_ResetParams = Academy.Instance.EnvironmentParameters;
    }

    /// <summary>
    /// Loop over body parts and reset them to initial conditions.
    /// </summary>
    public override void OnEpisodeBegin()
    {
        // Reset all of the body parts
        foreach (var bodyPart in m_JdController.bodyPartsDict.Values)
        {
            bodyPart.Reset(bodyPart);
        }

        // Random start rotation to help generalize
        hips.rotation = Quaternion.Euler(0, Random.Range(0.0f, 360.0f), 0);

        // Set our goal walking speed
        MTargetWalkingSpeed =
            randomizeWalkSpeedEachEpisode ? Random.Range(0.1f, m_maxWalkingSpeed) : MTargetWalkingSpeed;

        // Reset foot contact status tracking
        lastFootWasLeft = false;
        firstFootTouch = false;
    }


    /// <summary>
    /// Add relevant information on each body part to observations.
    /// </summary>
    public void CollectObservationBodyPart(BodyPart bp, VectorSensor sensor)
    {
        if (bp.rb.transform == footL || bp.rb.transform == footR)
        {
            bool touchingGround = bp.groundContact.touchingGround;
            sensor.AddObservation(touchingGround);

            if (touchingGround)
            {
                // Track the order of foot contact
                if (bp.rb.transform == footL)
                {
                    if (firstFootTouch)
                    {
                        if (lastFootWasLeft)
                        {
                            AddReward(-0.05f); // Penalize for not alternating
                        }
                        lastFootWasLeft = true;
                    }
                    else
                    {
                        firstFootTouch = true;
                        lastFootWasLeft = true;
                    }
                }
                else if (bp.rb.transform == footR)
                {
                    if (firstFootTouch)
                    {
                        if (!lastFootWasLeft)
                        {
                            AddReward(-0.05f); // Penalize for not alternating
                        }
                        lastFootWasLeft = false;
                    }
                    else
                    {
                        firstFootTouch = true;
                        lastFootWasLeft = false;
                    }
                }
            }
        }

        if (bp.rb.transform != footL && bp.rb.transform != footR && bp.rb.transform)
        {
            sensor.AddObservation(bp.rb.transform.localRotation);
        }
        if (bp.rb.transform != hips && bp.rb.transform)
        {
            sensor.AddObservation(bp.currentStrength / m_JdController.maxJointForceLimit);
        }
    }


    /// <summary>
    /// Loop over body parts to add them to observation.
    /// </summary>
    public override void CollectObservations(VectorSensor sensor)
    {
        foreach (var bodyPart in m_JdController.bodyPartsList)
        {
            CollectObservationBodyPart(bodyPart, sensor);
        }
    }

    public override void OnActionReceived(ActionBuffers actionBuffers)

    {
        var bpDict = m_JdController.bodyPartsDict;
        var i = -1;

        var continuousActions = actionBuffers.ContinuousActions;

        bpDict[thighL].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], 0);
        bpDict[thighR].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], 0);
        bpDict[shinL].SetJointTargetRotation(continuousActions[++i], 0, 0);
        bpDict[shinR].SetJointTargetRotation(continuousActions[++i], 0, 0);
        bpDict[footR].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], continuousActions[++i]);
        bpDict[footL].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], continuousActions[++i]);

        //update joint strength settings
        bpDict[thighL].SetJointStrength(continuousActions[++i]);
        bpDict[shinL].SetJointStrength(continuousActions[++i]);
        bpDict[footL].SetJointStrength(continuousActions[++i]);
        bpDict[thighR].SetJointStrength(continuousActions[++i]);
        bpDict[shinR].SetJointStrength(continuousActions[++i]);
        bpDict[footR].SetJointStrength(continuousActions[++i]);
    }

    void FixedUpdate()
    {
        CheckHipsHeight();
        // CheckHipsRotation();
        CheckDistance();
        LookAtTarget();
    }

    private void CheckHipsHeight()
    {
        float hipsHeight = hips.position.y;

        if (hipsHeight >= -0.3f) AddReward(+0.01f);
        else AddReward(-0.01f);
    }

    private void CheckHipsRotation()
    {
        float hipsRotation = hips.rotation.x;

        if (hipsRotation < 1.5 || hipsRotation > 2.5) AddReward(-0.01f);
        else AddReward(+0.03f);
    }

    private void CheckDistance()
    {
        float distanceToTarget = Vector3.Distance(hips.position, targetTransform.position);

        AddReward(distanceToTarget * 0.0001f);
    }

    private void LookAtTarget()
    {
        // Define the forward direction from the agent's hips
        Vector3 forwardDirection = hips.forward;

        // Cast a ray straight ahead from the agent
        Ray ray = new Ray(hips.position, forwardDirection);
        RaycastHit hit;

        // Check if the ray hits anything within range
        if (Physics.Raycast(ray, out hit, rayRange))
        {
            if (hit.transform == targetTransform)
            {
                AddReward(+1f);
                // Debug.Log("Target is directly ahead.");
            }
        }
    }
}