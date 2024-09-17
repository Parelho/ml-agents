using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SpaceStateLegs : MonoBehaviour
{
    int pernaAtual = 0;
    int seloPerna1 = 0;
    int seloPerna2 = 0;

    float g = 9.81f;
    float pi = Mathf.PI;

    float L1 = 0.0911f;
    float L3 = 0.3641f;
    float L4 = 0.3269f;
    float M1 = 4.55236f;
    float M3 = 1.86031f;
    float M4 = 0.25018f;
    float MY1 = 0;
    float MY3 = 0;
    float MY4 = 0;
    float MZ1;
    float MZ3;
    float MZ4;
    float XX1 = 0.15774f;
    float XX3 = 0.01156f;
    float XX4 = 0.00222f;

    float rVPP = 0.1f;
    float gam;

    float alpha0;
    float L0;
    float k = 10000;
    float c = 10;
    float offset = -0.04f;

    float[] qf = new float[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    float q31L = 0;
    float q32L = 0;
    float q41L = 0;
    float q42L = 0;
    float q1L = 0;
    float dq31L = 0;
    float dq32L = 0;
    float dq41L = 0;
    float dq42L = 0;
    float dq1L = 0;

    float[] L = new float[] { 0, 0, 0, 0, 0 };

    float pH1 = 0;
    float pH2 = 0;
    float pT1 = 0;
    float pT2 = 0;
    float pCoM1 = 0;
    float pCoM2 = 0;
    float pVPP1 = 0;
    float pVPP2 = 0;
    float pHead1 = 0;
    float pHead2 = 0;
    float pFem11 = 0;
    float pFem12 = 0;
    float pFem21 = 0;
    float pFem22 = 0;
    float pG11 = 0;
    float pG12 = 0;
    float pG21 = 0;
    float pG22 = 0;
    float pTib11 = 0;
    float pTib12 = 0;
    float pTib21 = 0;
    float pTib22 = 0;
    float pFoot11 = 0;
    float pFoot12 = 0;
    float pFoot21 = 0;
    float pFoot22 = 0;

    float[] xout = new float[14];

    void Start()
    {
        MZ1 = L1 / 2;
        MZ3 = L3 / 2;
        MZ4 = L4 / 2;
        gam = pi / 3;

        alpha0 = (pi * 115) / 180;
        L0 = L3 * Mathf.Sin((pi * 60) / 180) + L4 * Mathf.Sin((pi * 60) / 180);
    }

    float[] executaPerna(float q1L, float q31L, float q32L, float q41L, float q42L, float foot1, float foot2, float delta)
    {
        dq31L = (q31L - L[1]) / (delta);
        dq32L = (q32L - L[2]) / (delta);
        dq41L = (q41L - L[3]) / (delta);
        dq42L = (q42L - L[4]) / (delta);
        dq1L = (q1L - L[0]) / (delta);

        L = new float[] { q1L, q31L, q32L, q41L, q42L };

        if (pernaAtual == 0)
        {
            qf = new float[] { q32L, q31L, q42L, q41L, q1L, dq32L, dq31L, dq42L, dq41L, dq1L };
            xout = red2fullCoM5DoF(qf);
        }
        else
        {
            qf = new float[] { q31L, q32L, q41L, q42L, q1L, dq31L, dq32L, dq41L, dq42L, dq1L };
            xout = red2fullCoM5DoF(qf);
        }

        if (foot1 == 1 && foot2 == 0)
        {
            if (pernaAtual == 1 && seloPerna1 == 0)
            {
                pernaAtual = 0;
                seloPerna1 = 1;
                xout = impact5Links5DoFVPP(xout);
            }
        }
    }

    float[] red2fullCoM5DoF(float[] x)
    {
        float cmh = (M4 * MY4 * Mathf.Cos(x[0] + x[2] + x[4]) + M4 * MY4 * Mathf.Cos(x[1] + x[3] + x[4]) + L4 * M1 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[1] + x[3] + x[4]) + M3 * MY3 * Mathf.Cos(x[0] + x[4]) + M3 * MY3 * Mathf.Cos(x[1] + x[4]) + L3 * M1 * Mathf.Sin(x[0] + x[4]) + 2 * L3 * M3 * Mathf.Sin(x[0] + x[4]) + L3 * M4 * Mathf.Sin(x[0] + x[4]) - L3 * M4 * Mathf.Sin(x[1] + x[4]) - M3 * MZ3 * Mathf.Sin(x[0] + x[4]) - M3 * MZ3 * Mathf.Sin(x[1] + x[4]) + M1 * MY1 * Mathf.Cos(x[4]) - M1 * MZ1 * Mathf.Sin(x[4])) / (M1 + 2 * M3 + 2 * M4);
        float cmv = (M4 * MZ4 * Mathf.Cos(x[0] + x[2] + x[4]) - 2 * L4 * M3 * Mathf.Cos(x[0] + x[2] + x[4]) - 2 * L4 * M4 * Mathf.Cos(x[0] + x[2] + x[4]) - L4 * M1 * Mathf.Cos(x[0] + x[2] + x[4]) + M4 * MZ4 * Mathf.Cos(x[1] + x[3] + x[4]) + M4 * MY4 * Mathf.Sin(x[0] + x[2] + x[4]) + M4 * MY4 * Mathf.Sin(x[1] + x[3] + x[4]) - L3 * M1 * Mathf.Cos(x[0] + x[4]) - 2 * L3 * M3 * Mathf.Cos(x[0] + x[4]) - L3 * M4 * Mathf.Cos(x[0] + x[4]) + L3 * M4 * Mathf.Cos(x[1] + x[4]) + M3 * MZ3 * Mathf.Cos(x[0] + x[4]) + M3 * MZ3 * Mathf.Cos(x[1] + x[4]) + M3 * MY3 * Mathf.Sin(x[0] + x[4]) + M3 * MY3 * Mathf.Sin(x[1] + x[4]) + M1 * MZ1 * Mathf.Cos(x[4]) + M1 * MY1 * Mathf.Sin(x[4])) / (M1 + 2 * M3 + 2 * M4);
        float dcmh = (x[5] * (L4 * M1 * Mathf.Cos(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Cos(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Cos(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Cos(x[0] + x[2] + x[4]) - M4 * MY4 * Mathf.Sin(x[0] + x[2] + x[4]) + L3 * M1 * Mathf.Cos(x[0] + x[4]) + 2 * L3 * M3 * Mathf.Cos(x[0] + x[4]) + L3 * M4 * Mathf.Cos(x[0] + x[4]) - M3 * MZ3 * Mathf.Cos(x[0] + x[4]) - M3 * MY3 * Mathf.Sin(x[0] + x[4]))) / (M1 + 2 * M3 + 2 * M4) - (x[6] * (M4 * MZ4 * Mathf.Cos(x[1] + x[3] + x[4]) + M4 * MY4 * Mathf.Sin(x[1] + x[3] + x[4]) + L3 * M4 * Mathf.Cos(x[1] + x[4]) + M3 * MZ3 * Mathf.Cos(x[1] + x[4]) + M3 * MY3 * Mathf.Sin(x[1] + x[4]))) / (M1 + 2 * M3 + 2 * M4) - (x[8] * (M4 * MZ4 * Mathf.Cos(x[1] + x[3] + x[4]) + M4 * MY4 * Mathf.Sin(x[1] + x[3] + x[4]))) / (M1 + 2 * M3 + 2 * M4) - (x[9] * (M4 * MZ4 * Mathf.Cos(x[0] + x[2] + x[4]) - 2 * L4 * M3 * Mathf.Cos(x[0] + x[2] + x[4]) - 2 * L4 * M4 * Mathf.Cos(x[0] + x[2] + x[4]) - L4 * M1 * Mathf.Cos(x[0] + x[2] + x[4]) + M4 * MZ4 * Mathf.Cos(x[1] + x[3] + x[4]) + M4 * MY4 * Mathf.Sin(x[0] + x[2] + x[4]) + M4 * MY4 * Mathf.Sin(x[1] + x[3] + x[4]) - L3 * M1 * Mathf.Cos(x[0] + x[4]) - 2 * L3 * M3 * Mathf.Cos(x[0] + x[4]) - L3 * M4 * Mathf.Cos(x[0] + x[4]) + L3 * M4 * Mathf.Cos(x[1] + x[4]) + M3 * MZ3 * Mathf.Cos(x[0] + x[4]) + M3 * MZ3 * Mathf.Cos(x[1] + x[4]) + M3 * MY3 * Mathf.Sin(x[0] + x[4]) + M3 * MY3 * Mathf.Sin(x[1] + x[4]) + M1 * MZ1 * Mathf.Cos(x[4]) + M1 * MY1 * Mathf.Sin(x[4]))) / (M1 + 2 * M3 + 2 * M4) + (x[7] * (L4 * M1 * Mathf.Cos(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Cos(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Cos(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Cos(x[0] + x[2] + x[4]) - M4 * MY4 * Mathf.Sin(x[0] + x[2] + x[4]))) / (M1 + 2 * M3 + 2 * M4);
        float dcmv = (x[7] * (M4 * MY4 * Mathf.Cos(x[0] + x[2] + x[4]) + L4 * M1 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[0] + x[2] + x[4]))) / (M1 + 2 * M3 + 2 * M4) + (x[9] * (M4 * MY4 * Mathf.Cos(x[0] + x[2] + x[4]) + M4 * MY4 * Mathf.Cos(x[1] + x[3] + x[4]) + L4 * M1 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[1] + x[3] + x[4]) + M3 * MY3 * Mathf.Cos(x[0] + x[4]) + M3 * MY3 * Mathf.Cos(x[1] + x[4]) + L3 * M1 * Mathf.Sin(x[0] + x[4]) + 2 * L3 * M3 * Mathf.Sin(x[0] + x[4]) + L3 * M4 * Mathf.Sin(x[0] + x[4]) - L3 * M4 * Mathf.Sin(x[1] + x[4]) - M3 * MZ3 * Mathf.Sin(x[0] + x[4]) - M3 * MZ3 * Mathf.Sin(x[1] + x[4]) + M1 * MY1 * Mathf.Cos(x[4]) - M1 * MZ1 * Mathf.Sin(x[4]))) / (M1 + 2 * M3 + 2 * M4) + (x[5] * (M4 * MY4 * Mathf.Cos(x[0] + x[2] + x[4]) + L4 * M1 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M3 * Mathf.Sin(x[0] + x[2] + x[4]) + 2 * L4 * M4 * Mathf.Sin(x[0] + x[2] + x[4]) - M4 * MZ4 * Mathf.Sin(x[0] + x[2] + x[4]) + M3 * MY3 * Mathf.Cos(x[0] + x[4]) + L3 * M1 * Mathf.Sin(x[0] + x[4]) + 2 * L3 * M3 * Mathf.Sin(x[0] + x[4]) + L3 * M4 * Mathf.Sin(x[0] + x[4]) - M3 * MZ3 * Mathf.Sin(x[0] + x[4]))) / (M1 + 2 * M3 + 2 * M4) - (x[6] * (M4 * MZ4 * Mathf.Sin(x[1] + x[3] + x[4]) - M4 * MY4 * Mathf.Cos(x[1] + x[3] + x[4]) - M3 * MY3 * Mathf.Cos(x[1] + x[4]) + L3 * M4 * Mathf.Sin(x[1] + x[4]) + M3 * MZ3 * Mathf.Sin(x[1] + x[4]))) / (M1 + 2 * M3 + 2 * M4) + (x[8] * (M4 * MY4 * Mathf.Cos(x[1] + x[3] + x[4]) - M4 * MZ4 * Mathf.Sin(x[1] + x[3] + x[4]))) / (M1 + 2 * M3 + 2 * M4);

        float[] xCom = new float[] { x[0], x[1], x[2], x[3], x[4], cmh, cmv, x[5], x[6], x[7], x[8], x[9], dcmh, dcmv };

        return xCom;
    }

    float[] impact5Links5DoFVPP(float[] x)
    {
        float[] q = new float[] { x[0], x[1], x[2], x[3], x[4], x[5] };
        float[] qe = new float[] { x[0], x[1], x[2], x[3], x[4], x[5], x[6] };
        float[] dqeMinus = new float[] { x[7], x[8], x[9], x[10], x[11], x[12], x[13] };

        (float[,] De, float[,] E1e) = dynMod5Links7DoF(qe, dqeMinus);

        // Preenche E com as duas ultimas colunas de E1e
        float[,] E = new float[7, 2];
        int aux;
        for (int i = 0; i < E1e.GetLength(0); i++)
        {
            aux = 0;
            for (int j = 2; j < E1e.GetLength(1); j++)
            {
                E[i, aux++] = E1e[i, j];
            }
        }

        // Concatena os valores de De com -E horizontalmente
        float[,] M1 = new float[7, 9];
        for (int i = 0; i < 7; i++)
        {
            for (int j = 0; j < 7; j++)
            {
                M1[i, j] = De[i, j];
            }
            for (int j = 7; j < 9; j++)
            {
                M1[i, j] = -E[i, j - 7];
            }
        }
    }

    (float[,] D, float[,] E1) dynMod5Links7DoF(float[] qe, float[] dqeMinus)
    {
        float q31L = qe[0];
        float q32L = qe[1];
        float q41L = qe[2];
        float q42L = qe[3];
        float q1L = qe[4];

        float[,] D = new float[7, 7];
        D[0, 0] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) - 2 * L3 * M3 * M4 * MZ3 + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[0, 1] = -(L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) + M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) + M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) + M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) + 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[0, 2] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 + L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[0, 3] = -(M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[0, 4] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 - L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) - M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) - M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q31L + q41L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q31L + q41L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q31L + q41L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q31L + q41L) - L3 * M1 * M4 * MZ1 * Mathf.Cos(q31L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M1 * M3 * MY1 * MY3 * Mathf.Cos(q31L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q31L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY1 * Mathf.Sin(q31L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M1 * M3 * MY1 * MZ3 * Mathf.Sin(q31L) - M1 * M3 * MY3 * MZ1 * Mathf.Sin(q31L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[1, 0] = -(L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) + M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) + M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) + M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) + 3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) + 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[1, 1] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) - 2 * L3 * M3 * M4 * MZ3 + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[1, 2] = -(M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[1, 3] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 + L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[1, 4] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 - L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) - M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) - M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q32L + q42L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q32L + q42L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q32L + q42L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q32L + q42L) - L3 * M1 * M4 * MZ1 * Mathf.Cos(q32L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M1 * M3 * MY1 * MY3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY1 * Mathf.Sin(q32L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M1 * M3 * MY1 * MZ3 * Mathf.Sin(q32L) - M1 * M3 * MY3 * MZ1 * Mathf.Sin(q32L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[2, 0] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 + L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[2, 1] = -(M4 * (M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) + M3 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + L3 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) - M3 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) + M3 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L))) / (M1 + 2 * M3 + 2 * M4);
        D[2, 2] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4) / (M1 + 2 * M3 + 2 * M4);
        D[2, 3] = -(M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[2, 4] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q31L + q41L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q31L + q41L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q31L + q41L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q31L + q41L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[3, 0] = -(M4 * (M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) + M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + M3 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) + M3 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) - L3 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) + M3 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) - M3 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L))) / (M1 + 2 * M3 + 2 * M4);
        D[3, 1] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 + L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[3, 2] = -(M4 * M4 * Mathf.Cos(q31L - q32L + q41L - q42L) * (MY4 * MY4 + MZ4 * MZ4)) / (M1 + 2 * M3 + 2 * M4);
        D[3, 3] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4) / (M1 + 2 * M3 + 2 * M4);
        D[3, 4] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q32L + q42L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q32L + q42L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q32L + q42L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q32L + q42L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[4, 0] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 - L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) - M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) - M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q31L + q41L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q31L + q41L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q31L + q41L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q31L + q41L) - L3 * M1 * M4 * MZ1 * Mathf.Cos(q31L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M1 * M3 * MY1 * MY3 * Mathf.Cos(q31L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q31L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY1 * Mathf.Sin(q31L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M1 * M3 * MY1 * MZ3 * Mathf.Sin(q31L) - M1 * M3 * MY3 * MZ1 * Mathf.Sin(q31L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[4, 1] = (M1 * XX3 + M1 * XX4 + 2 * M3 * XX3 + 2 * M3 * XX4 + 2 * M4 * XX3 + 2 * M4 * XX4 + L3 * L3 * M4 * M4 + M3 * M3 * MY3 * MY3 + M4 * M4 * MY4 * MY4 + M3 * M3 * MZ3 * MZ3 + M4 * M4 * MZ4 * MZ4 + L3 * L3 * M1 * M4 + 2 * L3 * L3 * M3 * M4 + M1 * M3 * MY3 * MY3 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY3 * MY3 + 2 * M3 * M4 * MY4 * MY4 + M1 * M3 * MZ3 * MZ3 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ3 * MZ3 + 2 * M3 * M4 * MZ4 * MZ4 - L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) - M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) - M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q32L + q42L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q32L + q42L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q32L + q42L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q32L + q42L) - L3 * M1 * M4 * MZ1 * Mathf.Cos(q32L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M1 * M3 * MY1 * MY3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY1 * Mathf.Sin(q32L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M1 * M3 * MY1 * MZ3 * Mathf.Sin(q32L) - M1 * M3 * MY3 * MZ1 * Mathf.Sin(q32L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) - 2 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[4, 2] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q31L + q41L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q31L + q41L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q31L + q41L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q31L + q41L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L)) / (M1 + 2 * M3 + 2 * M4);
        D[4, 3] = (M1 * XX4 + 2 * M3 * XX4 + 2 * M4 * XX4 + M4 * M4 * MY4 * MY4 + M4 * M4 * MZ4 * MZ4 + M1 * M4 * MY4 * MY4 + 2 * M3 * M4 * MY4 * MY4 + M1 * M4 * MZ4 * MZ4 + 2 * M3 * M4 * MZ4 * MZ4 - M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) + L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) - L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - M1 * M4 * MY1 * MY4 * Mathf.Cos(q32L + q42L) - M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q32L + q42L) + M1 * M4 * MY1 * MZ4 * Mathf.Sin(q32L + q42L) - M1 * M4 * MY4 * MZ1 * Mathf.Sin(q32L + q42L) + L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 2 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) - M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 2 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L) - M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L)) / (M1 + 2 * M3 + 2 * M4);
        D[4, 4] = (M1 * XX1 + 2 * M1 * XX3 + 2 * M3 * XX1 + 2 * M1 * XX4 + 2 * M4 * XX1 + 4 * M3 * XX3 + 4 * M3 * XX4 + 4 * M4 * XX3 + 4 * M4 * XX4 + 2 * L3 * L3 * M4 * M4 + 2 * M3 * M3 * MY3 * MY3 + 2 * M4 * M4 * MY4 * MY4 + 2 * M3 * M3 * MZ3 * MZ3 + 2 * M4 * M4 * MZ4 * MZ4 + 2 * L3 * L3 * M1 * M4 + 4 * L3 * L3 * M3 * M4 + 2 * M1 * M3 * MY1 * MY1 + 2 * M1 * M4 * MY1 * MY1 + 2 * M1 * M3 * MY3 * MY3 + 2 * M1 * M4 * MY4 * MY4 + 4 * M3 * M4 * MY3 * MY3 + 4 * M3 * M4 * MY4 * MY4 + 2 * M1 * M3 * MZ1 * MZ1 + 2 * M1 * M4 * MZ1 * MZ1 + 2 * M1 * M3 * MZ3 * MZ3 + 2 * M1 * M4 * MZ4 * MZ4 + 4 * M3 * M4 * MZ3 * MZ3 + 4 * M3 * M4 * MZ4 * MZ4 - 2 * L3 * L3 * M4 * M4 * Mathf.Cos(q31L - q32L) - 2 * M3 * M3 * MY3 * MY3 * Mathf.Cos(q31L - q32L) - 2 * M3 * M3 * MZ3 * MZ3 * Mathf.Cos(q31L - q32L) - 2 * M4 * M4 * MY4 * MY4 * Mathf.Cos(q31L - q32L + q41L - q42L) - 2 * M4 * M4 * MZ4 * MZ4 * Mathf.Cos(q31L - q32L + q41L - q42L) - 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L + q41L) - 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L + q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q42L) - 2 * L3 * M4 * M4 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q42L) + 2 * L3 * M4 * M4 * MY4 * Mathf.Sin(q31L - q32L - q42L) - 4 * L3 * M3 * M4 * MZ3 - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L + q41L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L + q41L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L + q41L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L + q41L) - 2 * M1 * M4 * MY1 * MY4 * Mathf.Cos(q31L + q41L) - 2 * M1 * M4 * MY1 * MY4 * Mathf.Cos(q32L + q42L) - 2 * M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q31L + q41L) - 2 * M1 * M4 * MZ1 * MZ4 * Mathf.Cos(q32L + q42L) + 2 * M1 * M4 * MY1 * MZ4 * Mathf.Sin(q31L + q41L) - 2 * M1 * M4 * MY4 * MZ1 * Mathf.Sin(q31L + q41L) + 2 * M1 * M4 * MY1 * MZ4 * Mathf.Sin(q32L + q42L) - 2 * M1 * M4 * MY4 * MZ1 * Mathf.Sin(q32L + q42L) - 2 * L3 * M1 * M4 * MZ1 * Mathf.Cos(q31L) - 2 * L3 * M1 * M4 * MZ1 * Mathf.Cos(q32L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q41L) + 2 * L3 * M1 * M4 * MZ4 * Mathf.Cos(q42L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q41L) + 4 * L3 * M3 * M4 * MZ4 * Mathf.Cos(q42L) - 2 * M1 * M3 * MY1 * MY3 * Mathf.Cos(q31L) - 2 * M1 * M3 * MY1 * MY3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q41L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q42L) - 2 * M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q31L) - 2 * M1 * M3 * MZ1 * MZ3 * Mathf.Cos(q32L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q41L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q42L) - 2 * M3 * M4 * MY3 * MY4 * Mathf.Cos(q31L - q32L - q42L) - 2 * M3 * M4 * MZ3 * MZ4 * Mathf.Cos(q31L - q32L - q42L) + 2 * L3 * M1 * M4 * MY1 * Mathf.Sin(q31L) + 2 * L3 * M1 * M4 * MY1 * Mathf.Sin(q32L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q41L) + 2 * L3 * M1 * M4 * MY4 * Mathf.Sin(q42L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q41L) + 4 * L3 * M3 * M4 * MY4 * Mathf.Sin(q42L) + 2 * M1 * M3 * MY1 * MZ3 * Mathf.Sin(q31L) - 2 * M1 * M3 * MY3 * MZ1 * Mathf.Sin(q31L) + 2 * M1 * M3 * MY1 * MZ3 * Mathf.Sin(q32L) - 2 * M1 * M3 * MY3 * MZ1 * Mathf.Sin(q32L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q41L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q41L) + 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q42L) - 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q42L) - 2 * M3 * M4 * MY3 * MZ4 * Mathf.Sin(q31L - q32L - q42L) + 2 * M3 * M4 * MY4 * MZ3 * Mathf.Sin(q31L - q32L - q42L) - 4 * L3 * M3 * M4 * MZ3 * Mathf.Cos(q31L - q32L)) / (M1 + 2 * M3 + 2 * M4);
        D[5, 5] = M1 + 2 * M3 + 2 * M4;
        D[6, 6] = M1 + 2 * M3 + 2 * M4;

        float[,] E1 = new float[7, 4];
        E1[0, 0] = (M3 * (MZ3 * Mathf.Cos(q31L + q1L) + MY3 * Mathf.Sin(q31L + q1L)) + M4 * (L3 * Mathf.Cos(q31L + q1L) + MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Cos(q31L + q41L + q1L) - L3 * Mathf.Cos(q31L + q1L);
        E1[0, 1] = (M4 * (L3 * Mathf.Sin(q31L + q1L) - MY4 * Mathf.Cos(q31L + q41L + q1L) + MZ4 * Mathf.Sin(q31L + q41L + q1L)) - M3 * (MY3 * Mathf.Cos(q31L + q1L) - MZ3 * Mathf.Sin(q31L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Sin(q31L + q41L + q1L) - L3 * Mathf.Sin(q31L + q1L);
        E1[0, 2] = (M3 * (MZ3 * Mathf.Cos(q31L + q1L) + MY3 * Mathf.Sin(q31L + q1L)) + M4 * (L3 * Mathf.Cos(q31L + q1L) + MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[0, 3] = (M4 * (L3 * Mathf.Sin(q31L + q1L) - MY4 * Mathf.Cos(q31L + q41L + q1L) + MZ4 * Mathf.Sin(q31L + q41L + q1L)) - M3 * (MY3 * Mathf.Cos(q31L + q1L) - MZ3 * Mathf.Sin(q31L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[1, 0] = (M3 * (MZ3 * Mathf.Cos(q32L + q1L) + MY3 * Mathf.Sin(q32L + q1L)) + M4 * (L3 * Mathf.Cos(q32L + q1L) + MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[1, 1] = (M4 * (L3 * Mathf.Sin(q32L + q1L) - MY4 * Mathf.Cos(q32L + q42L + q1L) + MZ4 * Mathf.Sin(q32L + q42L + q1L)) - M3 * (MY3 * Mathf.Cos(q32L + q1L) - MZ3 * Mathf.Sin(q32L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[1, 2] = (M3 * (MZ3 * Mathf.Cos(q32L + q1L) + MY3 * Mathf.Sin(q32L + q1L)) + M4 * (L3 * Mathf.Cos(q32L + q1L) + MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Cos(q32L + q42L + q1L) - L3 * Mathf.Cos(q32L + q1L);
        E1[1, 3] = (M4 * (L3 * Mathf.Sin(q32L + q1L) - MY4 * Mathf.Cos(q32L + q42L + q1L) + MZ4 * Mathf.Sin(q32L + q42L + q1L)) - M3 * (MY3 * Mathf.Cos(q32L + q1L) - MZ3 * Mathf.Sin(q32L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Sin(q32L + q42L + q1L) - L3 * Mathf.Sin(q32L + q1L);
        E1[2, 0] = (M4 * (MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Cos(q31L + q41L + q1L);
        E1[2, 1] = -L4 * Mathf.Sin(q31L + q41L + q1L) - (M4 * (MY4 * Mathf.Cos(q31L + q41L + q1L) - MZ4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[2, 2] = (M4 * (MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[2, 3] = -(M4 * (MY4 * Mathf.Cos(q31L + q41L + q1L) - MZ4 * Mathf.Sin(q31L + q41L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[3, 0] = (M4 * (MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[3, 1] = -(M4 * (MY4 * Mathf.Cos(q32L + q42L + q1L) - MZ4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[3, 2] = (M4 * (MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Cos(q32L + q42L + q1L);
        E1[3, 3] = -L4 * Mathf.Sin(q32L + q42L + q1L) - (M4 * (MY4 * Mathf.Cos(q32L + q42L + q1L) - MZ4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4);
        E1[4, 0] = (M1 * (MZ1 * Mathf.Cos(q1L) + MY1 * Mathf.Sin(q1L)) + M3 * (MZ3 * Mathf.Cos(q31L + q1L) + MY3 * Mathf.Sin(q31L + q1L)) + M3 * (MZ3 * Mathf.Cos(q32L + q1L) + MY3 * Mathf.Sin(q32L + q1L)) + M4 * (L3 * Mathf.Cos(q31L + q1L) + MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L)) + M4 * (L3 * Mathf.Cos(q32L + q1L) + MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L3 * Mathf.Cos(q31L + q1L) - L4 * Mathf.Cos(q31L + q41L + q1L);
        E1[4, 1] = -L3 * Mathf.Sin(q31L + q1L) - (M1 * (MY1 * Mathf.Cos(q1L) - MZ1 * Mathf.Sin(q1L)) - M4 * (L3 * Mathf.Sin(q32L + q1L) - MY4 * Mathf.Cos(q32L + q42L + q1L) + MZ4 * Mathf.Sin(q32L + q42L + q1L)) - M4 * (L3 * Mathf.Sin(q31L + q1L) - MY4 * Mathf.Cos(q31L + q41L + q1L) + MZ4 * Mathf.Sin(q31L + q41L + q1L)) + M3 * (MY3 * Mathf.Cos(q31L + q1L) - MZ3 * Mathf.Sin(q31L + q1L)) + M3 * (MY3 * Mathf.Cos(q32L + q1L) - MZ3 * Mathf.Sin(q32L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Sin(q31L + q41L + q1L);
        E1[4, 2] = (M1 * (MZ1 * Mathf.Cos(q1L) + MY1 * Mathf.Sin(q1L)) + M3 * (MZ3 * Mathf.Cos(q31L + q1L) + MY3 * Mathf.Sin(q31L + q1L)) + M3 * (MZ3 * Mathf.Cos(q32L + q1L) + MY3 * Mathf.Sin(q32L + q1L)) + M4 * (L3 * Mathf.Cos(q31L + q1L) + MZ4 * Mathf.Cos(q31L + q41L + q1L) + MY4 * Mathf.Sin(q31L + q41L + q1L)) + M4 * (L3 * Mathf.Cos(q32L + q1L) + MZ4 * Mathf.Cos(q32L + q42L + q1L) + MY4 * Mathf.Sin(q32L + q42L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L3 * Mathf.Cos(q32L + q1L) - L4 * Mathf.Cos(q32L + q42L + q1L);
        E1[4, 3] = -L3 * Mathf.Sin(q32L + q1L) - (M1 * (MY1 * Mathf.Cos(q1L) - MZ1 * Mathf.Sin(q1L)) - M4 * (L3 * Mathf.Sin(q32L + q1L) - MY4 * Mathf.Cos(q32L + q42L + q1L) + MZ4 * Mathf.Sin(q32L + q42L + q1L)) - M4 * (L3 * Mathf.Sin(q31L + q1L) - MY4 * Mathf.Cos(q31L + q41L + q1L) + MZ4 * Mathf.Sin(q31L + q41L + q1L)) + M3 * (MY3 * Mathf.Cos(q31L + q1L) - MZ3 * Mathf.Sin(q31L + q1L)) + M3 * (MY3 * Mathf.Cos(q32L + q1L) - MZ3 * Mathf.Sin(q32L + q1L))) / (M1 + 2 * M3 + 2 * M4) - L4 * Mathf.Sin(q32L + q42L + q1L);
        E1[5, 0] = 1;
        E1[5, 2] = 1;
        E1[6, 1] = 1;
        E1[6, 3] = 1;

        return (D, E1);
    }
}
