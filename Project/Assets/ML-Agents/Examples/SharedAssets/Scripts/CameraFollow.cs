using UnityEngine;

namespace Unity.MLAgentsExamples
{
    public class CameraFollow : MonoBehaviour
    {
        [Tooltip("The target to follow")] 
        public Transform target;

        [Tooltip("The time it takes to move to the new position")]
        public float smoothingTime;

        [Tooltip("Optional additional offset for the camera (e.g., on z-axis)")]
        public Vector3 additionalOffset = new Vector3(0, 0, -11);

        private Vector3 m_Offset;
        private Vector3 m_CamVelocity;

        // Use this for initialization
        void Start()
        {
            // Calculate the initial offset
            m_Offset = gameObject.transform.position - target.position;
            
            // Apply additional offset if needed
            m_Offset += additionalOffset;
        }

        void FixedUpdate()
        {
            var newPosition = new Vector3(target.position.x + m_Offset.x, transform.position.y,
                target.position.z + m_Offset.z);

            gameObject.transform.position =
                Vector3.SmoothDamp(transform.position, newPosition, ref m_CamVelocity, smoothingTime, Mathf.Infinity,
                    Time.fixedDeltaTime);
        }
    }
}
