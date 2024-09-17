using System.Diagnostics;
using System.IO;
using UnityEngine;

public class RunPython : MonoBehaviour
{
    public string pythonPath = "/bin/python3.12";
    public string scriptPath = "executa_perna.py";
    private string filePath = "result.json";

    void Update()
    {
        RunPythonScript();

        string jsonContent = File.ReadAllText(filePath);
        ResultData result = JsonUtility.FromJson<ResultData>(jsonContent);
        UnityEngine.Debug.Log("t0: " + result.t0);
        UnityEngine.Debug.Log("t1: " + result.t1);
        UnityEngine.Debug.Log("t2: " + result.t2);
        UnityEngine.Debug.Log("t3: " + result.t3);
    }

    void RunPythonScript()
    {
        // Setup the process to run Python
        ProcessStartInfo psi = new ProcessStartInfo
        {
            FileName = pythonPath,
            Arguments = scriptPath,
            UseShellExecute = false,
            RedirectStandardError = true,
            CreateNoWindow = true
        };

        // Start the process
        Process process = Process.Start(psi);
        process.WaitForExit(); // Freeze frame until python is done executing
    }


    // Class to deserialize the JSON data
    [System.Serializable]
    public class ResultData
    {
        public float t0;
        public float t1;
        public float t2;
        public float t3;
    }
}
