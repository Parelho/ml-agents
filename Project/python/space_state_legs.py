# Versão mais recente da perna

import sys
import time
import math
import numpy as np


global perna_atual
global selo_perna_1 
global selo_perna_2
global torqueMax
global dcMax

perna_atual = 0
selo_perna_1 = 0
selo_perna_2 = 0


g = 9.81
pi = math.pi

L1 = 0.0911
L3 = 0.3641
L4 = 0.3269
M1 = 4.55236
M3 = 1.86031
M4 = 0.25018
MY1 = 0 #
MY3 = 0 #
MY4 = 0
MZ1 = L1/2
MZ3 = L3/2
MZ4 = L4/2
XX1 = 0.15774
XX3 = 0.01156
XX4 = 0.00222

rVPP = 0.1
gam = pi/3

alpha0 = (pi*115)/180 # 2.0071
L0 = L3*np.sin((pi*60)/180) + L4*np.sin((pi*60)/180) #  0.459
k = 10000
c = 10
offset = -0.04

q_f = np.array([0,0,0,0,0,0,0,0,0,0])

q31L=0
q32L=0
q41L=0
q42L=0
q1L= 0
dq31L=0
dq32L=0
dq41L=0
dq42L=0
dq1L= 0

L = [0,0,0,0,0]     


class Struct():
    pH1 = []
    pH2 = []
    pT1 = []
    pT2 = []
    pCoM1 = []
    pCoM2 = []
    pVPP1 = []
    pVPP2 = []
    pHead1 = []
    pHead2 = []
    pFem11 = []
    pFem12 = []
    pFem21 = []
    pFem22 = []
    pG11 = []
    pG12 = []
    pG21 = []
    pG22 = []
    pTib11 = []
    pTib12 = []
    pTib21 = []
    pTib22 = []
    pFoot11 = []
    pFoot12 = []
    pFoot21 = []
    pFoot22 = []

def jacobian_stance(q,dq):

    q32L=q[1]
    q41L=q[2]
    q42L=q[3]
    q1L=q[4]

    ST_length_T_S = (abs(L4*np.sin(q41L))**2 + abs(L3 + L4*np.cos(q41L))**2)**(1/2)
    ST_theta = q31L - pi/2 + np.arctan2(-(L3*np.cos(q31L + q1L)**2 + L3*np.sin(q31L + q1L)**2 + L4*np.cos(q31L + q41L + q1L)*np.cos(q31L + q1L) + L4*np.sin(q31L + q41L + q1L)*np.sin(q31L + q1L))/(np.cos(q31L + q1L)**2 + np.sin(q31L + q1L)**2), -(L4*np.cos(q31L + q41L + q1L)*np.sin(q31L + q1L) - L4*np.sin(q31L + q41L + q1L)*np.cos(q31L + q1L))/(np.cos(q31L + q1L)**2 + np.sin(q31L + q1L)**2))

    ST_J_polar_T_S = np.array(np.zeros((2,7)))
    ST_J_polar_T_S[0,2]=(2*L4*np.sign(L4*np.sin(q41L))*np.cos(q41L)*abs(L4*np.sin(q41L)) - 2*L4*abs(L3 + L4*np.cos(q41L))*np.sign(L3 + L4*np.cos(q41L))*np.sin(q41L))/(2*(abs(L4*np.sin(q41L))**2 + abs(L3 + L4*np.cos(q41L))**2)**(1/2))
    ST_J_polar_T_S[1,0] = 1
    ST_J_polar_T_S[1,2]=(2*np.imag(L4*np.cos(q41L))*np.real(L4*np.sin(q41L)) - 2*np.real(L4*np.cos(q41L))*np.imag(L4*np.sin(q41L)) + np.imag(L4*np.cos(q41L))**2 + np.real(L4*np.cos(q41L))**2 + np.imag(L4*np.sin(q41L))**2 + np.real(L4*np.sin(q41L))**2 + np.imag(L3)*np.imag(L4*np.cos(q41L)) + np.real(L3)*np.real(L4*np.cos(q41L)) + np.imag(L3)*np.real(L4*np.sin(q41L)) - np.real(L3)*np.imag(L4*np.sin(q41L)))/((np.imag(L3) + np.imag(L4*np.cos(q41L)))**2 - 2*np.imag(L4*np.sin(q41L))*(np.real(L3) + np.real(L4*np.cos(q41L))) + (np.real(L3) + np.real(L4*np.cos(q41L)))**2 + np.imag(L4*np.sin(q41L))**2 + np.real(L4*np.sin(q41L))**2 + 2*np.real(L4*np.sin(q41L))*(np.imag(L3) + np.imag(L4*np.cos(q41L))))

    SW_length_T_S = (abs(L4*np.sin(q42L))**2 + abs(L3 + L4*np.cos(q42L))**2)**(1/2)
    SW_theta = q32L - pi/2 + np.arctan2(-(L3*np.cos(q32L + q1L)**2 + L3*np.sin(q32L + q1L)**2 + L4*np.cos(q32L + q42L + q1L)*np.cos(q32L + q1L) + L4*np.sin(q32L + q42L + q1L)*np.sin(q32L + q1L))/(np.cos(q32L + q1L)**2 + np.sin(q32L + q1L)**2), -(L4*np.cos(q32L + q42L + q1L)*np.sin(q32L + q1L) - L4*np.sin(q32L + q42L + q1L)*np.cos(q32L + q1L))/(np.cos(q32L + q1L)**2 + np.sin(q32L + q1L)**2))

    SW_J_polar_T_S = np.array(np.zeros((2,7)))
    SW_J_polar_T_S[0,3]=(2*L4*np.sign(L4*np.sin(q42L))*np.cos(q42L)*abs(L4*np.sin(q42L)) - 2*L4*abs(L3 + L4*np.cos(q42L))*np.sign(L3 + L4*np.cos(q42L))*np.sin(q42L))/(2*(abs(L4*np.sin(q42L))**2 + abs(L3 + L4*np.cos(q42L))**2)**(1/2))
    SW_J_polar_T_S[1,1]=1
    SW_J_polar_T_S[1,3]=(2*np.imag(L4*np.cos(q42L))*np.real(L4*np.sin(q42L)) - 2*np.real(L4*np.cos(q42L))*np.imag(L4*np.sin(q42L)) + np.imag(L4*np.cos(q42L))**2 + np.real(L4*np.cos(q42L))**2 + np.imag(L4*np.sin(q42L))**2 + np.real(L4*np.sin(q42L))**2 + np.imag(L3)*np.imag(L4*np.cos(q42L)) + np.real(L3)*np.real(L4*np.cos(q42L)) + np.imag(L3)*np.real(L4*np.sin(q42L)) - np.real(L3)*np.imag(L4*np.sin(q42L)))/((np.imag(L3) + np.imag(L4*np.cos(q42L)))**2 - 2*np.imag(L4*np.sin(q42L))*(np.real(L3) + np.real(L4*np.cos(q42L))) + (np.real(L3) + np.real(L4*np.cos(q42L)))**2 + np.imag(L4*np.sin(q42L))**2 + np.real(L4*np.sin(q42L))**2 + 2*np.real(L4*np.sin(q42L))*(np.imag(L3) + np.imag(L4*np.cos(q42L))))

    return ST_length_T_S, ST_theta, ST_J_polar_T_S, SW_length_T_S, SW_theta,SW_J_polar_T_S

def angle2vec(vec1,vec2):
    x1 = vec1[0]
    y1 = vec1[1]
    x2 = vec2[0]
    y2 = vec2[1]
    angle = np.arctan2(x1*y2-y1*x2,x1*x2+y1*y2)
    return angle

def limb_position(x,pH_horiz):
    k = 5000
    out = Struct()
    nn = x.shape[0]
    if nn == 10:
        k == 0
    else:
        k=nn
    out.pH1 = pH_horiz
    out.pH2 = - L3*np.cos(x[0] + x[4]) - L4*np.cos(x[0] + x[2] + x[4])
    out.pT1 = out.pH1 + MY1*np.cos(x[4]) - MZ1*np.sin(x[4])
    out.pT2 = out.pH2 + MZ1*np.cos(x[4]) + MY1*np.sin(x[4])
    out.pCoM1 = out.pH1 + (M4*(L3*np.sin(x[0] + x[4]) - L3*np.sin(x[1] + x[4]) + MY4*np.cos(x[1] + x[3] + x[4]) + L4*np.sin(x[0] + x[2] + x[4]) - MZ4*np.sin(x[1] + x[3] + x[4])) + M1*(L3*np.sin(x[0] + x[4]) + MY1*np.cos(x[4]) - MZ1*np.sin(x[4]) + L4*np.sin(x[0] + x[2] + x[4])) + M4*(MY4*np.cos(x[0] + x[2] + x[4]) + L4*np.sin(x[0] + x[2] + x[4]) - MZ4*np.sin(x[0] + x[2] + x[4])) + M3*(MY3*np.cos(x[0] + x[4]) + L3*np.sin(x[0] + x[4]) - MZ3*np.sin(x[0] + x[4]) + L4*np.sin(x[0] + x[2] + x[4])) + M3*(MY3*np.cos(x[1] + x[4]) + L3*np.sin(x[0] + x[4]) - MZ3*np.sin(x[1] + x[4]) + L4*np.sin(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4) - L4*np.sin(x[0] + x[2] + x[4]) - L3*np.sin(x[0] + x[4])
    out.pCoM2 = out.pH2 + L3*np.cos(x[0] + x[4]) - (M1*(L3*np.cos(x[0] + x[4]) - MZ1*np.cos(x[4]) - MY1*np.sin(x[4]) + L4*np.cos(x[0] + x[2] + x[4])) - M4*(L3*np.cos(x[1] + x[4]) - L3*np.cos(x[0] + x[4]) - L4*np.cos(x[0] + x[2] + x[4]) + MZ4*np.cos(x[1] + x[3] + x[4]) + MY4*np.sin(x[1] + x[3] + x[4])) - M4*(MZ4*np.cos(x[0] + x[2] + x[4]) - L4*np.cos(x[0] + x[2] + x[4]) + MY4*np.sin(x[0] + x[2] + x[4])) + M3*(L3*np.cos(x[0] + x[4]) - MZ3*np.cos(x[0] + x[4]) - MY3*np.sin(x[0] + x[4]) + L4*np.cos(x[0] + x[2] + x[4])) + M3*(L3*np.cos(x[0] + x[4]) - MZ3*np.cos(x[1] + x[4]) - MY3*np.sin(x[1] + x[4]) + L4*np.cos(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4) + L4*np.cos(x[0] + x[2] + x[4])
    out.pVPP1 = out.pH1 + (M4*(L3*np.sin(x[0] + x[4]) - L3*np.sin(x[1] + x[4]) + MY4*np.cos(x[1] + x[3] + x[4]) + L4*np.sin(x[0] + x[2] + x[4]) - MZ4*np.sin(x[1] + x[3] + x[4])) + M1*(L3*np.sin(x[0] + x[4]) + MY1*np.cos(x[4]) - MZ1*np.sin(x[4]) + L4*np.sin(x[0] + x[2] + x[4])) + M4*(MY4*np.cos(x[0] + x[2] + x[4]) + L4*np.sin(x[0] + x[2] + x[4]) - MZ4*np.sin(x[0] + x[2] + x[4])) + M3*(MY3*np.cos(x[0] + x[4]) + L3*np.sin(x[0] + x[4]) - MZ3*np.sin(x[0] + x[4]) + L4*np.sin(x[0] + x[2] + x[4])) + M3*(MY3*np.cos(x[1] + x[4]) + L3*np.sin(x[0] + x[4]) - MZ3*np.sin(x[1] + x[4]) + L4*np.sin(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4) - rVPP*np.sin(gam) - L4*np.sin(x[0] + x[2] + x[4]) - L3*np.sin(x[0] + x[4])
    out.pVPP2 = out.pH2 + L3*np.cos(x[0] + x[4]) - (M1*(L3*np.cos(x[0] + x[4]) - MZ1*np.cos(x[4]) - MY1*np.sin(x[4]) + L4*np.cos(x[0] + x[2] + x[4])) - M4*(L3*np.cos(x[1] + x[4]) - L3*np.cos(x[0] + x[4]) - L4*np.cos(x[0] + x[2] + x[4]) + MZ4*np.cos(x[1] + x[3] + x[4]) + MY4*np.sin(x[1] + x[3] + x[4])) - M4*(MZ4*np.cos(x[0] + x[2] + x[4]) - L4*np.cos(x[0] + x[2] + x[4]) + MY4*np.sin(x[0] + x[2] + x[4])) + M3*(L3*np.cos(x[0] + x[4]) - MZ3*np.cos(x[0] + x[4]) - MY3*np.sin(x[0] + x[4]) + L4*np.cos(x[0] + x[2] + x[4])) + M3*(L3*np.cos(x[0] + x[4]) - MZ3*np.cos(x[1] + x[4]) - MY3*np.sin(x[1] + x[4]) + L4*np.cos(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4) + rVPP*np.cos(gam) + L4*np.cos(x[0] + x[2] + x[4])
    out.pHead1 = out.pH1 + -L1*np.sin(x[4])
    out.pHead2 = out.pH2 + L1*np.cos(x[4])
    out.pFem11 = out.pH1 + MY3*np.cos(x[0] + x[4]) - MZ3*np.sin(x[0] + x[4])
    out.pFem12 = out.pH2 + MZ3*np.cos(x[0] + x[4]) + MY3*np.sin(x[0] + x[4])
    out.pFem21 = out.pH1 + MY3*np.cos(x[1] + x[4]) - MZ3*np.sin(x[1] + x[4])
    out.pFem22 = out.pH2 + MZ3*np.cos(x[1] + x[4]) + MY3*np.sin(x[1] + x[4])
    out.pG11 = out.pH1 + -L3*np.sin(x[0] + x[4])
    out.pG12 = out.pH2 + L3*np.cos(x[0] + x[4])
    out.pG21 = out.pH1 + -L3*np.sin(x[1] + x[4])
    out.pG22 = out.pH2 + L3*np.cos(x[1] + x[4])
    out.pTib11 = out.pH1 + MY4*np.cos(x[0] + x[2] + x[4]) - L3*np.sin(x[0] + x[4]) - MZ4*np.sin(x[0] + x[2] + x[4])
    out.pTib12 = out.pH2 + L3*np.cos(x[0] + x[4]) + MZ4*np.cos(x[0] + x[2] + x[4]) + MY4*np.sin(x[0] + x[2] + x[4])
    out.pTib21 = out.pH1 + MY4*np.cos(x[1] + x[3] + x[4]) - L3*np.sin(x[1] + x[4]) - MZ4*np.sin(x[1] + x[3] + x[4])
    out.pTib22 = out.pH2 + L3*np.cos(x[1] + x[4]) + MZ4*np.cos(x[1] + x[3] + x[4]) + MY4*np.sin(x[1] + x[3] + x[4])
    out.pFoot11 = out.pH1 + - L3*np.sin(x[0] + x[4]) - L4*np.sin(x[0] + x[2] + x[4])
    out.pFoot12 = out.pH2 + L3*np.cos(x[0] + x[4]) + L4*np.cos(x[0] + x[2] + x[4])
    out.pFoot21 = out.pH1 + - L3*np.sin(x[1] + x[4]) - L4*np.sin(x[1] + x[3] + x[4])
    out.pFoot22 = out.pH2 + L3*np.cos(x[1] + x[4]) + L4*np.cos(x[1] + x[3] + x[4])
    return out

def swing_foot_height(x):
    nn = x.shape
    if nn != 1:
        h = L3*np.cos(x[1] + x[4]) - L3*np.cos(x[0] + x[4]) - L4*np.cos(x[0] + x[2] + x[4]) + L4*np.cos(x[1] + x[3] + x[4])
    else:
        h = L3*np.cos(x[2] + x[4]) - L3*np.cos(x[1] + x[4]) - L4*np.cos(x[1] + x[3] + x[4]) + L4*np.cos(x[2] + x[4] + x[4])
    return h

def Controller_VPP_store(q_f,dq_f):
    q_f_store = q_f
    dq_f_store = dq_f
    q = q_f_store[0:5]
    dq = dq_f_store[0:5]
    dq_f = dq_f_store.T
    out = limb_position(q.T,0)
    
    (ST_length_T_S,ST_theta,ST_J_polar_T_S,SW_length_T_S,SW_theta,SW_J_polar_T_S) = jacobian_stance(q,dq)
    ST_dX = np.dot(ST_J_polar_T_S,dq_f.T)
    ST_dL = ST_dX[0]
    
    SW_L = SW_length_T_S
    SW_J = SW_J_polar_T_S
    SW_dX = np.dot(SW_J,dq_f.T)
    SW_dL = SW_dX[0]
    SW_dth = SW_dX[1]
    
    ST_Fs = k*(L0-ST_length_T_S) - k/100*ST_dL
    if ST_Fs < 0:
        ST_Fs = 0 - k/100*ST_dL
    ST_B_L_C = np.array([[np.subtract(out.pCoM1,out.pFoot11)],[np.subtract(out.pCoM2,out.pFoot12)]])
    ST_B_L_H = np.array([[np.subtract(out.pH1,out.pFoot11)],[np.subtract(out.pH2,out.pFoot12)]])
    ST_zeta = angle2vec(ST_B_L_C,ST_B_L_H)
    cur_phi = q[4]-offset
    ST_cur_beta = -c*cur_phi
    ST_beta = ST_zeta + ST_cur_beta
    ST_Ft = ST_Fs*np.tan(ST_beta)
    ST_Tau_ff = ST_length_T_S*ST_Ft
    ST_F = np.array([[ST_Fs],[-ST_Tau_ff]] ,dtype=object)
    
    phi = cur_phi + pi/2
    SW_alpha = phi + SW_theta
    
    if ((SW_alpha - pi/2) < (pi*10)/180):
        SW_L_d = L0*4/5
    else:
        SW_L_d = L0
    
    SW_Fs = k*(SW_L_d-SW_L) - 100*SW_dL
    SW_Tau_ff = 0
    SW_Tau_fb = -5*2*(alpha0 - SW_alpha) + 1*SW_dth
    SW_F = np.array([[SW_Fs],[-(SW_Tau_ff+SW_Tau_fb)]])
    
    ST_F = np.array([[float(ST_F[0])],[float(ST_F[1])]])
    SW_F = np.array([[float(SW_F[0])],[float(SW_F[1])]])
    ST_tor = np.dot(ST_J_polar_T_S.T,ST_F)
    SW_tor = np.dot(SW_J.T,SW_F)
    D_fric = -np.diag([0.1, 0.1, 0.1, 0.1])
    a=[]
    a=np.delete(dq,-1)
    a=np.array([[a[0]],[a[1]],[a[2]],[a[3]]])
    u = np.array([[float(ST_tor[0])],[float(SW_tor[1])],[float(ST_tor[2])],[float(SW_tor[3])]]) +  np.dot(D_fric,a)
    return u

def red2full_CoM_5DoF(x):
    cmh = (M4*MY4*np.cos(x[0] + x[2] + x[4]) + M4*MY4*np.cos(x[1] + x[3] + x[4]) + L4*M1*np.sin(x[0] + x[2] + x[4]) + 2*L4*M3*np.sin(x[0] + x[2] + x[4]) + 2*L4*M4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[1] + x[3] + x[4]) + M3*MY3*np.cos(x[0] + x[4]) + M3*MY3*np.cos(x[1] + x[4]) + L3*M1*np.sin(x[0] + x[4]) + 2*L3*M3*np.sin(x[0] + x[4]) + L3*M4*np.sin(x[0] + x[4]) - L3*M4*np.sin(x[1] + x[4]) - M3*MZ3*np.sin(x[0] + x[4]) - M3*MZ3*np.sin(x[1] + x[4]) + M1*MY1*np.cos(x[4]) - M1*MZ1*np.sin(x[4]))/(M1 + 2*M3 + 2*M4)
    cmv = (M4*MZ4*np.cos(x[0] + x[2] + x[4]) - 2*L4*M3*np.cos(x[0] + x[2] + x[4]) - 2*L4*M4*np.cos(x[0] + x[2] + x[4]) - L4*M1*np.cos(x[0] + x[2] + x[4]) + M4*MZ4*np.cos(x[1] + x[3] + x[4]) + M4*MY4*np.sin(x[0] + x[2] + x[4]) + M4*MY4*np.sin(x[1] + x[3] + x[4]) - L3*M1*np.cos(x[0] + x[4]) - 2*L3*M3*np.cos(x[0] + x[4]) - L3*M4*np.cos(x[0] + x[4]) + L3*M4*np.cos(x[1] + x[4]) + M3*MZ3*np.cos(x[0] + x[4]) + M3*MZ3*np.cos(x[1] + x[4]) + M3*MY3*np.sin(x[0] + x[4]) + M3*MY3*np.sin(x[1] + x[4]) + M1*MZ1*np.cos(x[4]) + M1*MY1*np.sin(x[4]))/(M1 + 2*M3 + 2*M4)
    dcmh = (x[5]*(L4*M1*np.cos(x[0] + x[2] + x[4]) + 2*L4*M3*np.cos(x[0] + x[2] + x[4]) + 2*L4*M4*np.cos(x[0] + x[2] + x[4]) - M4*MZ4*np.cos(x[0] + x[2] + x[4]) - M4*MY4*np.sin(x[0] + x[2] + x[4]) + L3*M1*np.cos(x[0] + x[4]) + 2*L3*M3*np.cos(x[0] + x[4]) + L3*M4*np.cos(x[0] + x[4]) - M3*MZ3*np.cos(x[0] + x[4]) - M3*MY3*np.sin(x[0] + x[4])))/(M1 + 2*M3 + 2*M4) - (x[6]*(M4*MZ4*np.cos(x[1] + x[3] + x[4]) + M4*MY4*np.sin(x[1] + x[3] + x[4]) + L3*M4*np.cos(x[1] + x[4]) + M3*MZ3*np.cos(x[1] + x[4]) + M3*MY3*np.sin(x[1] + x[4])))/(M1 + 2*M3 + 2*M4) - (x[8]*(M4*MZ4*np.cos(x[1] + x[3] + x[4]) + M4*MY4*np.sin(x[1] + x[3] + x[4])))/(M1 + 2*M3 + 2*M4) - (x[9]*(M4*MZ4*np.cos(x[0] + x[2] + x[4]) - 2*L4*M3*np.cos(x[0] + x[2] + x[4]) - 2*L4*M4*np.cos(x[0] + x[2] + x[4]) - L4*M1*np.cos(x[0] + x[2] + x[4]) + M4*MZ4*np.cos(x[1] + x[3] + x[4]) + M4*MY4*np.sin(x[0] + x[2] + x[4]) + M4*MY4*np.sin(x[1] + x[3] + x[4]) - L3*M1*np.cos(x[0] + x[4]) - 2*L3*M3*np.cos(x[0] + x[4]) - L3*M4*np.cos(x[0] + x[4]) + L3*M4*np.cos(x[1] + x[4]) + M3*MZ3*np.cos(x[0] + x[4]) + M3*MZ3*np.cos(x[1] + x[4]) + M3*MY3*np.sin(x[0] + x[4]) + M3*MY3*np.sin(x[1] + x[4]) + M1*MZ1*np.cos(x[4]) + M1*MY1*np.sin(x[4])))/(M1 + 2*M3 + 2*M4) + (x[7]*(L4*M1*np.cos(x[0] + x[2] + x[4]) + 2*L4*M3*np.cos(x[0] + x[2] + x[4]) + 2*L4*M4*np.cos(x[0] + x[2] + x[4]) - M4*MZ4*np.cos(x[0] + x[2] + x[4]) - M4*MY4*np.sin(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4)
    dcmv = (x[7]*(M4*MY4*np.cos(x[0] + x[2] + x[4]) + L4*M1*np.sin(x[0] + x[2] + x[4]) + 2*L4*M3*np.sin(x[0] + x[2] + x[4]) + 2*L4*M4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[0] + x[2] + x[4])))/(M1 + 2*M3 + 2*M4) + (x[9]*(M4*MY4*np.cos(x[0] + x[2] + x[4]) + M4*MY4*np.cos(x[1] + x[3] + x[4]) + L4*M1*np.sin(x[0] + x[2] + x[4]) + 2*L4*M3*np.sin(x[0] + x[2] + x[4]) + 2*L4*M4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[1] + x[3] + x[4]) + M3*MY3*np.cos(x[0] + x[4]) + M3*MY3*np.cos(x[1] + x[4]) + L3*M1*np.sin(x[0] + x[4]) + 2*L3*M3*np.sin(x[0] + x[4]) + L3*M4*np.sin(x[0] + x[4]) - L3*M4*np.sin(x[1] + x[4]) - M3*MZ3*np.sin(x[0] + x[4]) - M3*MZ3*np.sin(x[1] + x[4]) + M1*MY1*np.cos(x[4]) - M1*MZ1*np.sin(x[4])))/(M1 + 2*M3 + 2*M4) + (x[5]*(M4*MY4*np.cos(x[0] + x[2] + x[4]) + L4*M1*np.sin(x[0] + x[2] + x[4]) + 2*L4*M3*np.sin(x[0] + x[2] + x[4]) + 2*L4*M4*np.sin(x[0] + x[2] + x[4]) - M4*MZ4*np.sin(x[0] + x[2] + x[4]) + M3*MY3*np.cos(x[0] + x[4]) + L3*M1*np.sin(x[0] + x[4]) + 2*L3*M3*np.sin(x[0] + x[4]) + L3*M4*np.sin(x[0] + x[4]) - M3*MZ3*np.sin(x[0] + x[4])))/(M1 + 2*M3 + 2*M4) - (x[6]*(M4*MZ4*np.sin(x[1] + x[3] + x[4]) - M4*MY4*np.cos(x[1] + x[3] + x[4]) - M3*MY3*np.cos(x[1] + x[4]) + L3*M4*np.sin(x[1] + x[4]) + M3*MZ3*np.sin(x[1] + x[4])))/(M1 + 2*M3 + 2*M4) + (x[8]*(M4*MY4*np.cos(x[1] + x[3] + x[4]) - M4*MZ4*np.sin(x[1] + x[3] + x[4])))/(M1 + 2*M3 + 2*M4)
    
    x_com = np.array([[float(x[0])],[float(x[1])],[float(x[2])],[float(x[3])],[float(x[4])],[float(cmh)],[float(cmv)],[float(x[5])],[float(x[6])],[float(x[7])],[float(x[8])],[float(x[9])],[float(dcmh)],[float(dcmv)]])
    return x_com

def impact_5links_5DoF_VPP(x):
    q = x[0:5]
    x_f = np.transpose([x])
    qe = x_f[0:7]
    dqe_minus = x_f[7:14]
    
    (De,E1e) = dyn_mod_5links_7DoF(qe,dqe_minus)
    E = E1e[:,2:4]
    M1 = np.concatenate((De,-E),axis=1)
    M2 = np.concatenate((E.T,np.zeros((2,2))),axis=1)
    A = np.concatenate((M1,M2),axis=0)
    B = np.concatenate((np.dot(De,dqe_minus),np.zeros((2,1))),axis=0)
    temp = np.linalg.inv(np.matrix(A))*B

    q_plus_new = np.array([[float(q[1])],[float(q[0])],[float(q[3])],[float(q[2])],[float(q[4])]])
    dq_plus_new = np.array([[float(temp[1])],[float(temp[0])],[float(temp[3])],[float(temp[2])],[float(temp[4])]])
    temp_a = red2full_CoM_5DoF(np.concatenate((q_plus_new,dq_plus_new),axis=0))
    return temp_a

def dyn_mod_5links_7DoF(q,dq):
    q31L=q[0]
    q32L=q[1]
    q41L=q[2]
    q42L=q[3]
    q1L=q[4]
    
    D=np.zeros((7,7))
    D[0,0]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 + 2*L3*M4**2*MZ4*np.cos(q41L) + 2*L3*M4**2*MY4*np.sin(q41L) - 2*L3*M3*M4*MZ3 + 2*L3*M1*M4*MZ4*np.cos(q41L) + 4*L3*M3*M4*MZ4*np.cos(q41L) - 2*M3*M4*MY3*MY4*np.cos(q41L) - 2*M3*M4*MZ3*MZ4*np.cos(q41L) + 2*L3*M1*M4*MY4*np.sin(q41L) + 4*L3*M3*M4*MY4*np.sin(q41L) + 2*M3*M4*MY3*MZ4*np.sin(q41L) - 2*M3*M4*MY4*MZ3*np.sin(q41L))/(M1 + 2*M3 + 2*M4)
    D[0,1]=-(L3**2*M4**2*np.cos(q31L - q32L) + M3**2*MY3**2*np.cos(q31L - q32L) + M3**2*MZ3**2*np.cos(q31L - q32L) + M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) + L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) - L3*M4**2*MY4*np.sin(q31L - q32L - q42L) + M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) + M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) + M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) + M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) + 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)                                                                  
    D[0,2]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 + L3*M4**2*MZ4*np.cos(q41L) + L3*M4**2*MY4*np.sin(q41L) + L3*M1*M4*MZ4*np.cos(q41L) + 2*L3*M3*M4*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q41L) - M3*M4*MZ3*MZ4*np.cos(q41L) + L3*M1*M4*MY4*np.sin(q41L) + 2*L3*M3*M4*MY4*np.sin(q41L) + M3*M4*MY3*MZ4*np.sin(q41L) - M3*M4*MY4*MZ3*np.sin(q41L))/(M1 + 2*M3 + 2*M4)
    D[0,3]=-(M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) - L3*M4**2*MY4*np.sin(q31L - q32L - q42L) + M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) + M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L))/(M1 + 2*M3 + 2*M4)
    D[0,4]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 - L3**2*M4**2*np.cos(q31L - q32L) - M3**2*MY3**2*np.cos(q31L - q32L) - M3**2*MZ3**2*np.cos(q31L - q32L) - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + 2*L3*M4**2*MZ4*np.cos(q41L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M4**2*MY4*np.sin(q41L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3 - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q31L + q41L) - M1*M4*MZ1*MZ4*np.cos(q31L + q41L) + M1*M4*MY1*MZ4*np.sin(q31L + q41L) - M1*M4*MY4*MZ1*np.sin(q31L + q41L) - L3*M1*M4*MZ1*np.cos(q31L) + 2*L3*M1*M4*MZ4*np.cos(q41L) + 4*L3*M3*M4*MZ4*np.cos(q41L) - M1*M3*MY1*MY3*np.cos(q31L) - 2*M3*M4*MY3*MY4*np.cos(q41L) - M1*M3*MZ1*MZ3*np.cos(q31L) - 2*M3*M4*MZ3*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY1*np.sin(q31L) + 2*L3*M1*M4*MY4*np.sin(q41L) + 4*L3*M3*M4*MY4*np.sin(q41L) + M1*M3*MY1*MZ3*np.sin(q31L) - M1*M3*MY3*MZ1*np.sin(q31L) + 2*M3*M4*MY3*MZ4*np.sin(q41L) - 2*M3*M4*MY4*MZ3*np.sin(q41L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)
    D[1,0]=-(L3**2*M4**2*np.cos(q31L - q32L) + M3**2*MY3**2*np.cos(q31L - q32L) + M3**2*MZ3**2*np.cos(q31L - q32L) + M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) + L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) - L3*M4**2*MY4*np.sin(q31L - q32L - q42L) + M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) + 3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) + M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) + M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) + 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)                                                                
    D[1,1]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 + 2*L3*M4**2*MZ4*np.cos(q42L) + 2*L3*M4**2*MY4*np.sin(q42L) - 2*L3*M3*M4*MZ3 + 2*L3*M1*M4*MZ4*np.cos(q42L) + 4*L3*M3*M4*MZ4*np.cos(q42L) - 2*M3*M4*MY3*MY4*np.cos(q42L) - 2*M3*M4*MZ3*MZ4*np.cos(q42L) + 2*L3*M1*M4*MY4*np.sin(q42L) + 4*L3*M3*M4*MY4*np.sin(q42L) + 2*M3*M4*MY3*MZ4*np.sin(q42L) - 2*M3*M4*MY4*MZ3*np.sin(q42L))/(M1 + 2*M3 + 2*M4)                                                                
    D[1,2]=-(M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) + L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) + M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L))/(M1 + 2*M3 + 2*M4)
    D[1,3]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 + L3*M4**2*MZ4*np.cos(q42L) + L3*M4**2*MY4*np.sin(q42L) + L3*M1*M4*MZ4*np.cos(q42L) + 2*L3*M3*M4*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q42L) - M3*M4*MZ3*MZ4*np.cos(q42L) + L3*M1*M4*MY4*np.sin(q42L) + 2*L3*M3*M4*MY4*np.sin(q42L) + M3*M4*MY3*MZ4*np.sin(q42L) - M3*M4*MY4*MZ3*np.sin(q42L))/(M1 + 2*M3 + 2*M4)
    D[1,4]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 - L3**2*M4**2*np.cos(q31L - q32L) - M3**2*MY3**2*np.cos(q31L - q32L) - M3**2*MZ3**2*np.cos(q31L - q32L) - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + 2*L3*M4**2*MZ4*np.cos(q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M4**2*MY4*np.sin(q42L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3 - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q32L + q42L) - M1*M4*MZ1*MZ4*np.cos(q32L + q42L) + M1*M4*MY1*MZ4*np.sin(q32L + q42L) - M1*M4*MY4*MZ1*np.sin(q32L + q42L) - L3*M1*M4*MZ1*np.cos(q32L) + 2*L3*M1*M4*MZ4*np.cos(q42L) + 4*L3*M3*M4*MZ4*np.cos(q42L) - M1*M3*MY1*MY3*np.cos(q32L) - 2*M3*M4*MY3*MY4*np.cos(q42L) - M1*M3*MZ1*MZ3*np.cos(q32L) - 2*M3*M4*MZ3*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY1*np.sin(q32L) + 2*L3*M1*M4*MY4*np.sin(q42L) + 4*L3*M3*M4*MY4*np.sin(q42L) + M1*M3*MY1*MZ3*np.sin(q32L) - M1*M3*MY3*MZ1*np.sin(q32L) + 2*M3*M4*MY3*MZ4*np.sin(q42L) - 2*M3*M4*MY4*MZ3*np.sin(q42L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)
    D[2,0]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 + L3*M4**2*MZ4*np.cos(q41L) + L3*M4**2*MY4*np.sin(q41L) + L3*M1*M4*MZ4*np.cos(q41L) + 2*L3*M3*M4*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q41L) - M3*M4*MZ3*MZ4*np.cos(q41L) + L3*M1*M4*MY4*np.sin(q41L) + 2*L3*M3*M4*MY4*np.sin(q41L) + M3*M4*MY3*MZ4*np.sin(q41L) - M3*M4*MY4*MZ3*np.sin(q41L))/(M1 + 2*M3 + 2*M4)
    D[2,1]=-(M4*(M4*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4*MZ4*np.cos(q31L - q32L + q41L) + M3*MY3*MY4*np.cos(q31L - q32L + q41L) + M3*MZ3*MZ4*np.cos(q31L - q32L + q41L) + L3*M4*MY4*np.sin(q31L - q32L + q41L) - M3*MY3*MZ4*np.sin(q31L - q32L + q41L) + M3*MY4*MZ3*np.sin(q31L - q32L + q41L)))/(M1 + 2*M3 + 2*M4)
    D[2,2]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2)/(M1 + 2*M3 + 2*M4)                                                      
    D[2,3]=-(M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L))/(M1 + 2*M3 + 2*M4)
    D[2,4]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + L3*M4**2*MZ4*np.cos(q41L) + L3*M4**2*MY4*np.sin(q41L) - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q31L + q41L) - M1*M4*MZ1*MZ4*np.cos(q31L + q41L) + M1*M4*MY1*MZ4*np.sin(q31L + q41L) - M1*M4*MY4*MZ1*np.sin(q31L + q41L) + L3*M1*M4*MZ4*np.cos(q41L) + 2*L3*M3*M4*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q41L) - M3*M4*MZ3*MZ4*np.cos(q41L) + L3*M1*M4*MY4*np.sin(q41L) + 2*L3*M3*M4*MY4*np.sin(q41L) + M3*M4*MY3*MZ4*np.sin(q41L) - M3*M4*MY4*MZ3*np.sin(q41L))/(M1 + 2*M3 + 2*M4)
    D[3,0]=-(M4*(M4*MY4**2*np.cos(q31L - q32L + q41L - q42L) + M4*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4*MZ4*np.cos(q31L - q32L - q42L) + M3*MY3*MY4*np.cos(q31L - q32L - q42L) + M3*MZ3*MZ4*np.cos(q31L - q32L - q42L) - L3*M4*MY4*np.sin(q31L - q32L - q42L) + M3*MY3*MZ4*np.sin(q31L - q32L - q42L) - M3*MY4*MZ3*np.sin(q31L - q32L - q42L)))/(M1 + 2*M3 + 2*M4)
    D[3,1]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 + L3*M4**2*MZ4*np.cos(q42L) + L3*M4**2*MY4*np.sin(q42L) + L3*M1*M4*MZ4*np.cos(q42L) + 2*L3*M3*M4*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q42L) - M3*M4*MZ3*MZ4*np.cos(q42L) + L3*M1*M4*MY4*np.sin(q42L) + 2*L3*M3*M4*MY4*np.sin(q42L) + M3*M4*MY3*MZ4*np.sin(q42L) - M3*M4*MY4*MZ3*np.sin(q42L))/(M1 + 2*M3 + 2*M4)
    D[3,2]=-(M4**2*np.cos(q31L - q32L + q41L - q42L)*(MY4**2 + MZ4**2))/(M1 + 2*M3 + 2*M4)
    D[3,3]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2)/(M1 + 2*M3 + 2*M4)                                                       
    D[3,4]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + L3*M4**2*MY4*np.sin(q42L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - M1*M4*MY1*MY4*np.cos(q32L + q42L) - M1*M4*MZ1*MZ4*np.cos(q32L + q42L) + M1*M4*MY1*MZ4*np.sin(q32L + q42L) - M1*M4*MY4*MZ1*np.sin(q32L + q42L) + L3*M1*M4*MZ4*np.cos(q42L) + 2*L3*M3*M4*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q42L) - M3*M4*MZ3*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY4*np.sin(q42L) + 2*L3*M3*M4*MY4*np.sin(q42L) + M3*M4*MY3*MZ4*np.sin(q42L) - M3*M4*MY4*MZ3*np.sin(q42L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L))/(M1 + 2*M3 + 2*M4)
    D[4,0]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 - L3**2*M4**2*np.cos(q31L - q32L) - M3**2*MY3**2*np.cos(q31L - q32L) - M3**2*MZ3**2*np.cos(q31L - q32L) - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + 2*L3*M4**2*MZ4*np.cos(q41L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M4**2*MY4*np.sin(q41L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3 - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q31L + q41L) - M1*M4*MZ1*MZ4*np.cos(q31L + q41L) + M1*M4*MY1*MZ4*np.sin(q31L + q41L) - M1*M4*MY4*MZ1*np.sin(q31L + q41L) - L3*M1*M4*MZ1*np.cos(q31L) + 2*L3*M1*M4*MZ4*np.cos(q41L) + 4*L3*M3*M4*MZ4*np.cos(q41L) - M1*M3*MY1*MY3*np.cos(q31L) - 2*M3*M4*MY3*MY4*np.cos(q41L) - M1*M3*MZ1*MZ3*np.cos(q31L) - 2*M3*M4*MZ3*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY1*np.sin(q31L) + 2*L3*M1*M4*MY4*np.sin(q41L) + 4*L3*M3*M4*MY4*np.sin(q41L) + M1*M3*MY1*MZ3*np.sin(q31L) - M1*M3*MY3*MZ1*np.sin(q31L) + 2*M3*M4*MY3*MZ4*np.sin(q41L) - 2*M3*M4*MY4*MZ3*np.sin(q41L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)
    D[4,1]=(M1*XX3 + M1*XX4 + 2*M3*XX3 + 2*M3*XX4 + 2*M4*XX3 + 2*M4*XX4 + L3**2*M4**2 + M3**2*MY3**2 + M4**2*MY4**2 + M3**2*MZ3**2 + M4**2*MZ4**2 + L3**2*M1*M4 + 2*L3**2*M3*M4 + M1*M3*MY3**2 + M1*M4*MY4**2 + 2*M3*M4*MY3**2 + 2*M3*M4*MY4**2 + M1*M3*MZ3**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ3**2 + 2*M3*M4*MZ4**2 - L3**2*M4**2*np.cos(q31L - q32L) - M3**2*MY3**2*np.cos(q31L - q32L) - M3**2*MZ3**2*np.cos(q31L - q32L) - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + 2*L3*M4**2*MZ4*np.cos(q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M4**2*MY4*np.sin(q42L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3 - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q32L + q42L) - M1*M4*MZ1*MZ4*np.cos(q32L + q42L) + M1*M4*MY1*MZ4*np.sin(q32L + q42L) - M1*M4*MY4*MZ1*np.sin(q32L + q42L) - L3*M1*M4*MZ1*np.cos(q32L) + 2*L3*M1*M4*MZ4*np.cos(q42L) + 4*L3*M3*M4*MZ4*np.cos(q42L) - M1*M3*MY1*MY3*np.cos(q32L) - 2*M3*M4*MY3*MY4*np.cos(q42L) - M1*M3*MZ1*MZ3*np.cos(q32L) - 2*M3*M4*MZ3*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY1*np.sin(q32L) + 2*L3*M1*M4*MY4*np.sin(q42L) + 4*L3*M3*M4*MY4*np.sin(q42L) + M1*M3*MY1*MZ3*np.sin(q32L) - M1*M3*MY3*MZ1*np.sin(q32L) + 2*M3*M4*MY3*MZ4*np.sin(q42L) - 2*M3*M4*MY4*MZ3*np.sin(q42L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) - 2*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)
    D[4,2]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + L3*M4**2*MZ4*np.cos(q41L) + L3*M4**2*MY4*np.sin(q41L) - M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - M1*M4*MY1*MY4*np.cos(q31L + q41L) - M1*M4*MZ1*MZ4*np.cos(q31L + q41L) + M1*M4*MY1* MZ4*np.sin(q31L + q41L) - M1*M4*MY4*MZ1*np.sin(q31L + q41L) + L3*M1*M4*MZ4*np.cos(q41L) + 2*L3*M3*M4*MZ4*np.cos(q41L) - M3*M4*MY3*MY4*np.cos(q41L) - M3*M4*MZ3*MZ4*np.cos(q41L) + L3*M1*M4*MY4*np.sin(q41L) + 2*L3*M3*M4*MY4*np.sin(q41L) + M3*M4*MY3*MZ4*np.sin(q41L) - M3*M4*MY4*MZ3*np.sin(q41L))/(M1 + 2*M3 + 2*M4)
    D[4,3]=(M1*XX4 + 2*M3*XX4 + 2*M4*XX4 + M4**2*MY4**2 + M4**2*MZ4**2 + M1*M4*MY4**2 + 2*M3*M4*MY4**2 + M1*M4*MZ4**2 + 2*M3*M4*MZ4**2 - M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) + L3*M4**2*MZ4*np.cos(q42L) - L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + L3*M4**2*MY4*np.sin(q42L) + L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - M1*M4*MY1*MY4*np.cos(q32L + q42L) - M1*M4*MZ1*MZ4*np.cos(q32L + q42L) + M1*M4*MY1*MZ4*np.sin(q32L + q42L) - M1*M4*MY4*MZ1*np.sin(q32L + q42L) + L3*M1*M4*MZ4*np.cos(q42L) + 2*L3*M3*M4*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q42L) - M3*M4*MZ3*MZ4*np.cos(q42L) - M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + L3*M1*M4*MY4*np.sin(q42L) + 2*L3*M3*M4*MY4*np.sin(q42L) + M3*M4*MY3*MZ4*np.sin(q42L) - M3*M4*MY4*MZ3*np.sin(q42L) - M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L))/(M1 + 2*M3 + 2*M4)
    D[4,4]=(M1*XX1 + 2*M1*XX3 + 2*M3*XX1 + 2*M1*XX4 + 2*M4*XX1 + 4*M3*XX3 + 4*M3*XX4 + 4*M4*XX3 + 4*M4*XX4 + 2*L3**2*M4**2 + 2*M3**2*MY3**2 + 2*M4**2*MY4**2 + 2*M3**2*MZ3**2 + 2*M4**2*MZ4**2 + 2*L3**2*M1*M4 + 4*L3**2*M3*M4 + 2*M1*M3*MY1**2 + 2*M1*M4*MY1**2 + 2*M1*M3*MY3**2 + 2*M1*M4*MY4**2 + 4*M3*M4*MY3**2 + 4*M3*M4*MY4**2 + 2*M1*M3*MZ1**2 + 2*M1*M4*MZ1**2 + 2*M1*M3*MZ3**2 + 2*M1*M4*MZ4**2 + 4*M3*M4*MZ3**2 + 4*M3*M4*MZ4**2 - 2*L3**2*M4**2*np.cos(q31L - q32L) - 2*M3**2*MY3**2*np.cos(q31L - q32L) - 2*M3**2*MZ3**2*np.cos(q31L - q32L) - 2*M4**2*MY4**2*np.cos(q31L - q32L + q41L - q42L) - 2*M4**2*MZ4**2*np.cos(q31L - q32L + q41L - q42L) - 2*L3*M4**2*MZ4*np.cos(q31L - q32L + q41L) - 2*L3*M4**2*MY4*np.sin(q31L - q32L + q41L) + 2*L3*M4**2*MZ4*np.cos(q41L) + 2*L3*M4**2*MZ4*np.cos(q42L) - 2*L3*M4**2*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M4**2*MY4*np.sin(q41L) + 2*L3*M4**2*MY4*np.sin(q42L) + 2*L3*M4**2*MY4*np.sin(q31L - q32L - q42L) - 4*L3*M3*M4*MZ3 - 2*M3*M4*MY3*MY4*np.cos(q31L - q32L + q41L) - 2*M3*M4*MZ3*MZ4*np.cos(q31L - q32L + q41L) + 2*M3*M4*MY3*MZ4*np.sin(q31L - q32L + q41L) - 2*M3*M4*MY4*MZ3*np.sin(q31L - q32L + q41L) - 2*M1*M4*MY1*MY4*np.cos(q31L + q41L) - 2*M1*M4*MY1*MY4*np.cos(q32L + q42L) - 2*M1*M4*MZ1*MZ4*np.cos(q31L + q41L) - 2*M1*M4*MZ1*MZ4*np.cos(q32L + q42L) + 2*M1*M4*MY1*MZ4*np.sin(q31L + q41L) - 2*M1*M4*MY4*MZ1*np.sin(q31L + q41L) + 2*M1*M4*MY1*MZ4*np.sin(q32L + q42L) - 2*M1*M4*MY4*MZ1*np.sin(q32L + q42L) - 2*L3*M1*M4*MZ1*np.cos(q31L) - 2*L3*M1*M4*MZ1*np.cos(q32L) + 2*L3*M1*M4*MZ4*np.cos(q41L) + 2*L3*M1*M4*MZ4*np.cos(q42L) + 4*L3*M3*M4*MZ4*np.cos(q41L) + 4*L3*M3*M4*MZ4*np.cos(q42L) - 2*M1*M3*MY1*MY3*np.cos(q31L) - 2*M1*M3*MY1*MY3*np.cos(q32L) - 2*M3*M4*MY3*MY4*np.cos(q41L) - 2*M3*M4*MY3*MY4*np.cos(q42L) - 2*M1*M3*MZ1*MZ3*np.cos(q31L) - 2*M1*M3*MZ1*MZ3*np.cos(q32L) - 2*M3*M4*MZ3*MZ4*np.cos(q41L) - 2*M3*M4*MZ3*MZ4*np.cos(q42L) - 2*M3*M4*MY3*MY4*np.cos(q31L - q32L - q42L) - 2*M3*M4*MZ3*MZ4*np.cos(q31L - q32L - q42L) + 2*L3*M1*M4*MY1*np.sin(q31L) + 2*L3*M1*M4*MY1*np.sin(q32L) + 2*L3*M1*M4*MY4*np.sin(q41L) + 2*L3*M1*M4*MY4*np.sin(q42L) + 4*L3*M3*M4*MY4*np.sin(q41L) + 4*L3*M3*M4*MY4*np.sin(q42L) + 2*M1*M3*MY1*MZ3*np.sin(q31L) - 2*M1*M3*MY3*MZ1*np.sin(q31L) + 2*M1*M3*MY1*MZ3*np.sin(q32L) - 2*M1*M3*MY3*MZ1*np.sin(q32L) + 2*M3*M4*MY3*MZ4*np.sin(q41L) - 2*M3*M4*MY4*MZ3*np.sin(q41L) + 2*M3*M4*MY3*MZ4*np.sin(q42L) - 2*M3*M4*MY4*MZ3*np.sin(q42L) - 2*M3*M4*MY3*MZ4*np.sin(q31L - q32L - q42L) + 2*M3*M4*MY4*MZ3*np.sin(q31L - q32L - q42L) - 4*L3*M3*M4*MZ3*np.cos(q31L - q32L))/(M1 + 2*M3 + 2*M4)
    D[5,5]=M1 + 2*M3 + 2*M4
    D[6,6]=M1 + 2*M3 + 2*M4
    
    E1=np.zeros((7,4))
    E1[0,0]=(M3*(MZ3*np.cos(q31L + q1L) + MY3*np.sin(q31L + q1L)) + M4*(L3*np.cos(q31L + q1L) + MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.cos(q31L + q41L + q1L) - L3*np.cos(q31L + q1L)
    E1[0,1]=(M4*(L3*np.sin(q31L + q1L) - MY4*np.cos(q31L + q41L + q1L) + MZ4*np.sin(q31L + q41L + q1L)) - M3*(MY3*np.cos(q31L + q1L) - MZ3*np.sin(q31L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.sin(q31L + q41L + q1L) - L3*np.sin(q31L + q1L)
    E1[0,2]=(M3*(MZ3*np.cos(q31L + q1L) + MY3*np.sin(q31L + q1L)) + M4*(L3*np.cos(q31L + q1L) + MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[0,3]=(M4*(L3*np.sin(q31L + q1L) - MY4*np.cos(q31L + q41L + q1L) + MZ4*np.sin(q31L + q41L + q1L)) - M3*(MY3*np.cos(q31L + q1L) - MZ3*np.sin(q31L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[1,0]=(M3*(MZ3*np.cos(q32L + q1L) + MY3*np.sin(q32L + q1L)) + M4*(L3*np.cos(q32L + q1L) + MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[1,1]=(M4*(L3*np.sin(q32L + q1L) - MY4*np.cos(q32L + q42L + q1L) + MZ4*np.sin(q32L + q42L + q1L)) - M3*(MY3*np.cos(q32L + q1L) - MZ3*np.sin(q32L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[1,2]=(M3*(MZ3*np.cos(q32L + q1L) + MY3*np.sin(q32L + q1L)) + M4*(L3*np.cos(q32L + q1L) + MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.cos(q32L + q42L + q1L) - L3*np.cos(q32L + q1L)
    E1[1,3]=(M4*(L3*np.sin(q32L + q1L) - MY4*np.cos(q32L + q42L + q1L) + MZ4*np.sin(q32L + q42L + q1L)) - M3*(MY3*np.cos(q32L + q1L) - MZ3*np.sin(q32L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.sin(q32L + q42L + q1L) - L3*np.sin(q32L + q1L)
    E1[2,0]=(M4*(MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.cos(q31L + q41L + q1L)
    E1[2,1]=-L4*np.sin(q31L + q41L + q1L) - (M4*(MY4*np.cos(q31L + q41L + q1L) - MZ4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[2,2]=(M4*(MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[2,3]=-(M4*(MY4*np.cos(q31L + q41L + q1L) - MZ4*np.sin(q31L + q41L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[3,0]=(M4*(MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[3,1]=-(M4*(MY4*np.cos(q32L + q42L + q1L) - MZ4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[3,2]=(M4*(MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.cos(q32L + q42L + q1L)
    E1[3,3]=-L4*np.sin(q32L + q42L + q1L) - (M4*(MY4*np.cos(q32L + q42L + q1L) - MZ4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4)
    E1[4,0]=(M1*(MZ1*np.cos(q1L) + MY1*np.sin(q1L)) + M3*(MZ3*np.cos(q31L + q1L) + MY3*np.sin(q31L + q1L)) + M3*(MZ3*np.cos(q32L + q1L) + MY3*np.sin(q32L + q1L)) + M4*(L3*np.cos(q31L + q1L) + MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)) + M4*(L3*np.cos(q32L + q1L) + MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4) - L3*np.cos(q31L + q1L) - L4*np.cos(q31L + q41L + q1L)
    E1[4,1]=-L3*np.sin(q31L + q1L) - (M1*(MY1*np.cos(q1L) - MZ1*np.sin(q1L)) - M4*(L3*np.sin(q32L + q1L) - MY4*np.cos(q32L + q42L + q1L) + MZ4*np.sin(q32L + q42L + q1L)) - M4*(L3*np.sin(q31L + q1L) - MY4*np.cos(q31L + q41L + q1L) + MZ4*np.sin(q31L + q41L + q1L)) + M3*(MY3*np.cos(q31L + q1L) - MZ3*np.sin(q31L + q1L)) + M3*(MY3*np.cos(q32L + q1L) - MZ3*np.sin(q32L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.sin(q31L + q41L + q1L)
    E1[4,2]=(M1*(MZ1*np.cos(q1L) + MY1*np.sin(q1L)) + M3*(MZ3*np.cos(q31L + q1L) + MY3*np.sin(q31L + q1L)) + M3*(MZ3*np.cos(q32L + q1L) + MY3*np.sin(q32L + q1L)) + M4*(L3*np.cos(q31L + q1L) + MZ4*np.cos(q31L + q41L + q1L) + MY4*np.sin(q31L + q41L + q1L)) + M4*(L3*np.cos(q32L + q1L) + MZ4*np.cos(q32L + q42L + q1L) + MY4*np.sin(q32L + q42L + q1L)))/(M1 + 2*M3 + 2*M4) - L3*np.cos(q32L + q1L) - L4*np.cos(q32L + q42L + q1L)
    E1[4,3]=-L3*np.sin(q32L + q1L) - (M1*(MY1*np.cos(q1L) - MZ1*np.sin(q1L)) - M4*(L3*np.sin(q32L + q1L) - MY4*np.cos(q32L + q42L + q1L) + MZ4*np.sin(q32L + q42L + q1L)) - M4*(L3*np.sin(q31L + q1L) - MY4*np.cos(q31L + q41L + q1L) + MZ4*np.sin(q31L + q41L + q1L)) + M3*(MY3*np.cos(q31L + q1L) - MZ3*np.sin(q31L + q1L)) + M3*(MY3*np.cos(q32L + q1L) - MZ3*np.sin(q32L + q1L)))/(M1 + 2*M3 + 2*M4) - L4*np.sin(q32L + q42L + q1L)
    E1[5,0]=1
    E1[5,2]=1
    E1[6,1]=1
    E1[6,3]=1
    return D, E1

def mapRange(valor): 
    torqueMax = 10
    dcMax = 20

    if valor > 0:
        if valor > torqueMax:
            valor = torqueMax        
        return (valor/torqueMax) * dcMax

    if valor <= 0:
        if valor < -torqueMax:
            valor = -torqueMax        
        return (valor/torqueMax) * dcMax


def executa_perna(q1L, q31L, q32L, q41L, q42L, foot1, foot2,delta): 
    
    global perna_atual
    global selo_perna_1 
    global selo_perna_2 
    global L
    # q1L                   #self.data5[0,-1]/ 1000   tronco
    # q31L                  #self.data1[0,-1]/ 1000   femur swing (esquerdo)
    # q32L                  #self.data2[0,-1]/ 1000 #femur stance (direito)
    # q41L                  #self.data3[0,-1]/ 1000 #tibia swing (esquerdo)
    # q42L                   #self.data4[0,-1]/ 1000 #tibia stance (direito)
    # foot1 
    # foot2 

    dq31L =  (q31L - L[1])/(delta) # DEBUBMARTINS
    dq32L =  (q32L - L[2])/(delta)
    dq41L =  (q41L - L[3])/(delta)
    dq42L =  (q42L - L[4])/(delta)
    dq1L  =   (q1L - L[0])/(delta)

    # dq31L =  (q31L - L[1])
    # dq32L =  (q32L - L[2])
    # dq41L =  (q41L - L[3])
    # dq42L =  (q42L - L[4])
    # dq1L  =   (q1L - L[0])

    L = [q1L,q31L,q32L,q41L,q42L]
    #print(selo_perna_1)
    #print(selo_perna_2)
    # print(f' p_atual: {perna_atual} ', end="")

    if (perna_atual == 0): #perna swing esquerda
        q_f = np.array([q32L,q31L,q42L,q41L,q1L,dq32L,dq31L,dq42L,dq41L,dq1L])
        x0 = red2full_CoM_5DoF(q_f)
        xout = x0.T
        xout_f = xout[0]
        x1 = impact_5links_5DoF_VPP(xout_f)
    else: #perna swing direita (perna_atual == 1)
        q_f = np.array([q31L,q32L,q41L,q42L,q1L,dq31L,dq32L,dq41L,dq42L,dq1L])        
        x0 = red2full_CoM_5DoF(q_f)
        xout = x0.T
        xout_f = xout[0]

    if (( foot1 == 1) & ( foot2 == 0 )): # limiar do sensor perna direita ativo em 1
        if ((perna_atual == 1) & (selo_perna_1 == 0)): #precisa trocar de perna, pois perna_atual=1 eh a perna2
            perna_atual = 0 #troca a perna swing da direita pra esquerda
            selo_perna_1 = 1 #trava a perna direita
            x1 = impact_5links_5DoF_VPP(xout_f)
            xout = x1.T
            xout_f = xout[0]

        else:
            selo_perna_1 = 0  #destrava acionamento pois a perna nao esta no chao

    if (( foot1 == 0) & ( foot2 == 1 )): #limiar do sensor de perna esquerda ativo em 1
        if ((perna_atual == 0) & (selo_perna_2 == 0)):  # precisa trocar de perna, pois perna_atual=0 eh a perna1
            perna_atual = 1 #troca a perna swing da esquerda pra direita
            selo_perna_2 = 1 #trava a perna esquerda
            x1 = impact_5links_5DoF_VPP(xout_f)
            xout = x1.T
            xout_f = xout[0]
        else:
            selo_perna_2 = 0

    torqueout = Controller_VPP_store(xout_f[0:7].T,xout_f[7:14].T)

    u1 =  float(torqueout[0])      #motor superior perna direita 
    u2 = -float(torqueout[1])      #motor inferior perna direita 
    u3 = -float(torqueout[2])      #motor superior perna esquerda 
    u4 =  float(torqueout[3])      #motor inferior perna esquerda
        
    return [mapRange(u1), mapRange(u2), mapRange(u3), mapRange(u4)]
