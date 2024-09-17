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
    [Range(2, 10)]
    [SerializeField]
    //The walking speed to try and achieve
    private float m_TargetWalkingSpeed = 10;

    public float MTargetWalkingSpeed // property
    {
        get { return m_TargetWalkingSpeed; }
        set { m_TargetWalkingSpeed = Mathf.Clamp(value, 2f, m_maxWalkingSpeed); }
    }

    const float m_maxWalkingSpeed = 10; //The max walking speed

    //Should the agent sample a new goal velocity each episode?
    //If true, walkSpeed will be randomly set between zero and m_maxWalkingSpeed in OnEpisodeBegin() 
    //If false, the goal velocity will be walkingSpeed
    public bool randomizeWalkSpeedEachEpisode;

    //The direction an agent will walk during training.
    private Vector3 m_WorldDirToWalk = Vector3.left;

    [Header("Target To Walk Towards")] public Transform target; //Target the agent will walk towards during training.

    [Header("Body Parts")] public Transform hips;
    public Transform thighL;
    public Transform shinL;
    public Transform footL;
    public Transform thighR;
    public Transform shinR;
    public Transform footR;

    private bool isLeftFootOnGround;
    private bool isRightFootOnGround;
    private bool lastFootWasLeft;
    private bool firstFootTouch;



    //This will be used as a stabilized model space reference point for observations
    //Because ragdolls can move erratically during training, using a stabilized reference transform improves learning
    OrientationCubeController m_OrientationCube;

    //The indicator graphic gameobject that points towards the target
    DirectionIndicator m_DirectionIndicator;
    JointDriveController m_JdController;
    EnvironmentParameters m_ResetParams;

    public override void Initialize()
    {
        m_OrientationCube = GetComponentInChildren<OrientationCubeController>();
        m_DirectionIndicator = GetComponentInChildren<DirectionIndicator>();

        //Setup each body part
        m_JdController = GetComponent<JointDriveController>();
        m_JdController.SetupBodyPart(hips);
        m_JdController.SetupBodyPart(thighL);
        m_JdController.SetupBodyPart(shinL);
        m_JdController.SetupBodyPart(footL);
        m_JdController.SetupBodyPart(thighR);
        m_JdController.SetupBodyPart(shinR);
        m_JdController.SetupBodyPart(footR);

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

        UpdateOrientationObjects();

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
        // GROUND CHECK
        if (bp.rb.transform == footL || bp.rb.transform == footR)
        {
            bool touchingGround = bp.groundContact.touchingGround;
            sensor.AddObservation(touchingGround);

            // // Update foot contact status
            if (bp.rb.transform == footL)
            {
                isLeftFootOnGround = touchingGround;
            }
            else if (bp.rb.transform == footR)
            {
                isRightFootOnGround = touchingGround;
            }

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

        // Get velocities in the context of our orientation cube's space
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformDirection(bp.rb.velocity));
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformDirection(bp.rb.angularVelocity));

        // Get position relative to hips in the context of our orientation cube's space
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformDirection(bp.rb.position - hips.position));

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
        var cubeForward = m_OrientationCube.transform.forward;

        //velocity we want to match
        var velGoal = cubeForward * MTargetWalkingSpeed;
        //ragdoll's avg vel
        var avgVel = GetAvgVelocity();

        //current ragdoll velocity. normalized 
        sensor.AddObservation(Vector3.Distance(velGoal, avgVel));
        //avg body vel relative to cube
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformDirection(avgVel));
        //vel goal relative to cube
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformDirection(velGoal));

        //rotation deltas
        sensor.AddObservation(Quaternion.FromToRotation(hips.forward, cubeForward));

        //Position of target position relative to cube
        sensor.AddObservation(m_OrientationCube.transform.InverseTransformPoint(target.transform.position));

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

    //Update OrientationCube and DirectionIndicator
    void UpdateOrientationObjects()
    {
        m_WorldDirToWalk = target.position - hips.position;
        m_OrientationCube.UpdateOrientation(hips, target);
        if (m_DirectionIndicator)
        {
            m_DirectionIndicator.MatchOrientation(m_OrientationCube.transform);
        }
    }

    void FixedUpdate()
    {
        UpdateOrientationObjects();

        var cubeForward = m_OrientationCube.transform.forward;

        // Set reward for this step according to mixture of the following elements.
        // a. Match target speed
        var matchSpeedReward = GetMatchingVelocityReward(cubeForward * MTargetWalkingSpeed, GetAvgVelocity());

        // Check for NaNs
        if (float.IsNaN(matchSpeedReward))
        {
            throw new ArgumentException(
                "NaN in moveTowardsTargetReward.\n" +
                $" cubeForward: {cubeForward}\n" +
                $" hips.velocity: {m_JdController.bodyPartsDict[hips].rb.velocity}\n" +
                $" maximumWalkingSpeed: {m_maxWalkingSpeed}"
            );
        }

        // b. Rotation alignment with target direction.
        var lookAtTargetReward = (Vector3.Dot(cubeForward, hips.forward) + 1) * .5F;

        // Check for NaNs
        if (float.IsNaN(lookAtTargetReward))
        {
            throw new ArgumentException(
                "NaN in lookAtTargetReward.\n" +
                $" cubeForward: {cubeForward}\n" +
                $" head.forward: {hips.forward}"
            );
        }

        // Apply rewards
        if (lookAtTargetReward >= 0.75)
        {
            AddReward(+lookAtTargetReward * matchSpeedReward);
        }
        else
        {
            AddReward(-lookAtTargetReward * matchSpeedReward);
        }
        if (matchSpeedReward > 0.9)
        {
            AddReward(+0.5f);
        }

        CheckHipsHeight();
    }

    private void CheckHipsHeight()
    {
        float hipsHeight = hips.position.y;

        if (hipsHeight >= -0.1f) AddReward(+0.03f);
        else AddReward(-0.03f);
    }

    //Returns the average velocity of all of the body parts
    //Using the velocity of the hips only has shown to result in more erratic movement from the limbs, so...
    //...using the average helps prevent this erratic movement
    Vector3 GetAvgVelocity()
    {
        Vector3 velSum = Vector3.zero;
        Vector3 avgVel = Vector3.zero;

        //ALL RBS
        int numOfRB = 0;
        foreach (var item in m_JdController.bodyPartsList)
        {
            numOfRB++;
            velSum += item.rb.velocity;
        }

        avgVel = velSum / numOfRB;
        return avgVel;
    }

    //normalized value of the difference in avg speed vs goal walking speed.
    public float GetMatchingVelocityReward(Vector3 velocityGoal, Vector3 actualVelocity)
    {
        //distance between our actual velocity and goal velocity
        var velDeltaMagnitude = Mathf.Clamp(Vector3.Distance(actualVelocity, velocityGoal), 0, MTargetWalkingSpeed);

        //return the value on a declining sigmoid shaped curve that decays from 1 to 0
        //This reward will approach 1 if it matches perfectly and approach zero as it deviates
        return Mathf.Pow(1 - Mathf.Pow(velDeltaMagnitude / MTargetWalkingSpeed, 2), 2);
    }

    /// <summary>
    /// Agent touched the target
    /// </summary>
    public void TouchedTarget()
    {
        AddReward(1f);
    }
}

// using System;
// using UnityEngine;
// using Unity.MLAgents;
// using Unity.MLAgents.Actuators;
// using Unity.MLAgentsExamples;
// using Unity.MLAgents.Sensors;
// using BodyPart = Unity.MLAgentsExamples.BodyPart;
// using Random = UnityEngine.Random;

// public class WalkerAgent : Agent
// {
//     [Header("Walk Speed")]
//     [Range(1, 10)]
//     [SerializeField]
//     private float m_TargetWalkingSpeed = 10;

//     public float MTargetWalkingSpeed
//     {
//         get { return m_TargetWalkingSpeed; }
//         set { m_TargetWalkingSpeed = Mathf.Clamp(value, 1f, m_maxWalkingSpeed); }
//     }

//     const float m_maxWalkingSpeed = 10;

//     public bool randomizeWalkSpeedEachEpisode;
//     [Header("Body Parts")] public Transform hips;
//     public Transform thighL;
//     public Transform shinL;
//     public Transform footL;
//     public Transform thighR;
//     public Transform shinR;
//     public Transform footR;
//     public Transform targetTransform;

//     private bool lastFootWasLeft;
//     private bool firstFootTouch;

//     private Rigidbody hipsRb;
//     private float nextPushTime = 0f;

//     JointDriveController m_JdController;
//     EnvironmentParameters m_ResetParams;

//     public override void Initialize()
//     {
//         m_JdController = GetComponent<JointDriveController>();
//         m_JdController.SetupBodyPart(hips);
//         m_JdController.SetupBodyPart(thighL);
//         m_JdController.SetupBodyPart(shinL);
//         m_JdController.SetupBodyPart(footL);
//         m_JdController.SetupBodyPart(thighR);
//         m_JdController.SetupBodyPart(shinR);
//         m_JdController.SetupBodyPart(footR);

//         targetTransform = transform.root.Find("Target");

//         hipsRb = hips.GetComponent<Rigidbody>();

//         m_ResetParams = Academy.Instance.EnvironmentParameters;
//     }

//     public override void OnEpisodeBegin()
//     {
//         foreach (var bodyPart in m_JdController.bodyPartsDict.Values)
//         {
//             bodyPart.Reset(bodyPart);
//         }

//         hips.rotation = Quaternion.Euler(0, 0, 0);

//         MTargetWalkingSpeed =
//             randomizeWalkSpeedEachEpisode ? Random.Range(0.1f, m_maxWalkingSpeed) : MTargetWalkingSpeed;

//         lastFootWasLeft = false;
//         firstFootTouch = false;
//     }

//     public void CollectObservationBodyPart(BodyPart bp, VectorSensor sensor)
//     {
//         if (bp.rb.transform == footL || bp.rb.transform == footR)
//         {
//             bool touchingGround = bp.groundContact.touchingGround;
//             sensor.AddObservation(touchingGround);

//             if (touchingGround)
//             {
//                 if (bp.rb.transform == footL)
//                 {
//                     if (firstFootTouch)
//                     {
//                         if (lastFootWasLeft)
//                         {
//                             AddReward(-0.01f);
//                         }
//                         lastFootWasLeft = true;
//                     }
//                     else
//                     {
//                         firstFootTouch = true;
//                         lastFootWasLeft = true;
//                     }
//                 }
//                 else if (bp.rb.transform == footR)
//                 {
//                     if (firstFootTouch)
//                     {
//                         if (!lastFootWasLeft)
//                         {
//                             AddReward(-0.01f);
//                         }
//                         lastFootWasLeft = false;
//                     }
//                     else
//                     {
//                         firstFootTouch = true;
//                         lastFootWasLeft = false;
//                     }
//                 }
//             }
//         }

//         if (bp.rb.transform != footL && bp.rb.transform != footR && bp.rb.transform)
//         {
//             sensor.AddObservation(bp.rb.transform.localRotation);
//         }
//         if (bp.rb.transform != hips && bp.rb.transform)
//         {
//             sensor.AddObservation(bp.currentStrength / m_JdController.maxJointForceLimit);
//         }
//     }

//     public override void CollectObservations(VectorSensor sensor)
//     {
//         foreach (var bodyPart in m_JdController.bodyPartsList)
//         {
//             CollectObservationBodyPart(bodyPart, sensor);
//         }

//         // Add the walking speed observation
//         sensor.AddObservation(MTargetWalkingSpeed);

//         // Calculate and add normalized distance to target
//         float maxDistance = 12f; // Assuming a reasonable max distance for normalization
//         float distanceToTarget = Vector3.Distance(hips.position, targetTransform.position);
//         float normalizedDistance = Mathf.Clamp01(distanceToTarget / maxDistance);
//         sensor.AddObservation(normalizedDistance);
//     }

//     public override void OnActionReceived(ActionBuffers actionBuffers)
//     {
//         var bpDict = m_JdController.bodyPartsDict;
//         var i = -1;

//         var continuousActions = actionBuffers.ContinuousActions;

//         bpDict[thighL].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], 0);
//         bpDict[thighR].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], 0);
//         bpDict[shinL].SetJointTargetRotation(continuousActions[++i], 0, 0);
//         bpDict[shinR].SetJointTargetRotation(continuousActions[++i], 0, 0);
//         bpDict[footR].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], continuousActions[++i]);
//         bpDict[footL].SetJointTargetRotation(continuousActions[++i], continuousActions[++i], continuousActions[++i]);

//         // update joint strength settings
//         bpDict[thighL].SetJointStrength(continuousActions[++i]);
//         bpDict[shinL].SetJointStrength(continuousActions[++i]);
//         bpDict[footL].SetJointStrength(continuousActions[++i]);
//         bpDict[thighR].SetJointStrength(continuousActions[++i]);
//         bpDict[shinR].SetJointStrength(continuousActions[++i]);
//         bpDict[footR].SetJointStrength(continuousActions[++i]);
//     }

//     void FixedUpdate()
//     {
//         CheckHipsHeight();
//         CheckDistance();
//         LookAtTarget();
//     }

//     private void CheckHipsHeight()
//     {
//         float hipsHeight = hips.position.y;

//         if (hipsHeight >= -0.1f) AddReward(+0.001f);
//         else AddReward(-0.001f);
//     }

//     private void CheckDistance()
//     {
//         float maxDistance = 12f;
//         float distanceToTarget = Vector3.Distance(hips.position, targetTransform.position);

//         float normalizedDistance = Mathf.Clamp01(1 - (distanceToTarget / maxDistance));

//         AddReward(normalizedDistance * normalizedDistance);
//     }

//     private void LookAtTarget()
//     {
//         float currentRotationY = hips.rotation.eulerAngles.y;

//         float rotationDifference = Mathf.DeltaAngle(currentRotationY, 0f);
//         float maxRotationDifference = 180f;
//         float normalizedRotation = Mathf.Clamp01(1 - (Mathf.Abs(rotationDifference) / maxRotationDifference));

//         AddReward(normalizedRotation * 0.001f);
//     }

//     int pushDirection;
//     Vector3 force;
//     private void PushHips()
//     {
//         if (Time.time >= nextPushTime)
//         {
//             pushDirection = Random.Range(0, 3);

//             force = pushDirection == 0 ? hips.forward : pushDirection == 1 ? hips.forward * -1 : pushDirection == 2 ? hips.right : hips.right * -1;
//             hipsRb.AddForce(force * 250);

//             nextPushTime = Time.time + Random.Range(5f, 10f);
//         }
//     }
// }
