#!/usr/bin/env python
# coding: utf-8

#Riyaz Ghulam Anwary
#12217052

#TM4112
#Project Reservoir Simulation 2020/2021
#Due date: 23/12/2020

# In[264]:

import pandas as pd
import numpy as np
from scipy.linalg import solve


# In[265]:


def readdata():
    global nx, ny, nz, tx, ty, tz, Pi, Swi
    global phi0, crock, p_ref, kx, ky, kz, nrock
    global Sw, Krw, Kro, Pcow, ros, rgs, rws, npvto, npvtg
    global Rs, Poil, Bo, Muo, Pgas, Bg, Mug, Pw_ref, Bw_ref, Cw, Muw_ref, Vscw 
    global Nw, wlx, wly, wlz , wrv, wr
    with open("Data.txt", "r") as rr:
        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=int)
        nx = temp[0]
        ny = temp[1]
        nz = temp[2]

        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        tx = temp[0]
        ty = temp[1]
        tz = temp[2]

        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        Pi = temp[0]
        Swi = temp[1]

        for i in range(0, 6):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        phi0 = temp[0]
        crock = temp[1]
        p_ref = temp[2]

        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        kx = temp[0]
        ky = temp[1]
        kz = temp[2]

        for i in range(0, 4):
            rr.readline()
        
        nrock = int(rr.readline())
        Sw = np.zeros(nrock, dtype=float)
        Krw = np.zeros(nrock, dtype=float)
        Kro = np.zeros(nrock, dtype=float)
        Pcow = np.zeros(nrock, dtype=float)

        for i in range(0, 6):
            rr.readline()

        for i in range(0, nrock):
            line = rr.readline()
            temp = np.array(line.split(), dtype=float)
            Sw[i] = temp[0]
            Krw[i] = temp[1]
            Kro[i] = temp[2]
            Pcow[i] = temp[3]
        
        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        ros = temp[0]
        rgs = temp[1]
        rws = temp[2]

        for i in range(0, 5):
            rr.readline()
        
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        npvto = int(temp[0])
        npvtg = int(temp[1])
        
        for i in range(0, 7):
            rr.readline()

        Rs = np.zeros(npvto, dtype=float)
        Poil = np.zeros(npvto, dtype=float)
        Bo = np.zeros(npvto, dtype=float)
        Muo = np.zeros(npvto, dtype=float)
        
        for i in range(0, npvto):
            line = rr.readline()
            temp = np.array(line.split(), dtype=float)
            Rs[i] = temp[0]
            Poil[i] = temp[1]
            Bo[i] = temp[2]
            Muo[i] = temp[3]
        
        for i in range(0, 6):
            rr.readline()

        Pgas = np.zeros(npvtg, dtype=float)
        Bg = np.zeros(npvtg, dtype=float)
        Mug = np.zeros(npvtg, dtype=float)
        
        for i in range(0, npvtg):
            line = rr.readline()
            temp = np.array(line.split(), dtype=float)
            Pgas[i] = temp[0]
            Bg[i] = temp[1]
            Mug[i] = temp[2]
        
        for i in range(0, 8):
            rr.readline()
         
        line = rr.readline()
        temp = np.array(line.split(), dtype=float)
        Pw_ref = temp[0]
        Bw_ref = temp[1]
        Cw = temp[2]
        Muw_ref = temp[3]
        Vscw = temp[4]
        
        for i in range(0, 4):
            rr.readline()
        
        Nw = int(rr.readline()) # numOfWell
        
        wlx = np.zeros(Nw, dtype=int)
        wly = np.zeros(Nw, dtype=int)
        wlz = np.zeros(Nw, dtype=int)

        for i in range(0, 4):
            rr.readline()
        
        for i in range(0, Nw):
            line = rr.readline()
            temp = np.array(line.split(), dtype=int)
            wlx[i] = temp[0]-1
            wly[i] = temp[1]-1
            wlz[i] = temp[2]-1
        
        for i in range(0, 4):
            rr.readline()
        
        wrv = np.zeros(Nw, dtype=int)
        for i in range(0, Nw):
            wrv[i] = int(rr.readline())
        
        for i in range(0, 4):
            rr.readline()
        
        wr = np.zeros((Nw, np.amax(wrv), 2), dtype=float)
        for i in range(0, Nw):
            rr.readline()
            for j in range(0, wrv[i]):
                line = rr.readline()
                temp = np.array(line.split(), dtype=float)
                wr[i][j][0] = temp[0] #time
                wr[i][j][1] = temp[1] #rate

        return


# In[266]:


def interpolate(tabx, taby, x):
    i = 0
    if (x<tabx[i]):
        y = taby[i]
    elif (x>tabx[-1]):
        y = taby[-1]
    else:
        x1 = tabx[i]
        x2 = tabx[i+1]
        y1 = taby[i]
        y2 = taby[i+1]
        while(x>tabx[i+1]):
            i += 1
            x1 = tabx[i]
            x2 = tabx[i + 1]
            y1 = taby[i]
            y2 = taby[i + 1]
        y = y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)
    return y


# In[267]:


#fungsi subroutine fluid property
def bo (p): #Fungsi Formation Volume Factor
    return interpolate (Poil, Bo, p)

def miuo (p):
    return interpolate (Poil, Muo, p)

def rs (p):
    return interpolate (Poil, Rs, p) * 1000/5.6146

def bw (p):
    x = Cw * (p-Pw_ref)
    y = Bw_ref / (1+x+(x**2)/2)
    return y

def miuw (p):
    Y = -1 * Cw * (p-Pw_ref)
    Miuw = Muw_ref / (1+Y+(Y**2)/2) 
    return Miuw

def bg (p):
    return interpolate(Pgas, Bg, p)*1000*5.6146

def rho_oil (p):
    rhoo = (ros + rs(p)*rgs) / bo(p) / 144
    return rhoo

def rho_gas (p):
    rhog = rgs / bg(p) / 144
    return rhog

def rho_water (p):
    rhow = rws / bw(p) / 144
    return rhow

#fungsi subroutine turunan fluid property
def dbop (p):
    z = 0.001
    v1 = bo(p)
    v2 = bo(p + z)
    return (v2 - v1)/z

def dmiuop (p):
    z = 0.001 #interval
    v1 = miuo(p)
    v2 = miuo(p + z)
    return (v2 - v1)/z 

def drsp (p):
    z = 0.001
    v1 = rs(p)
    v2 = rs(p + z)
    return (v2 - v1)/z

def dbwp (p):
    z = 0.001 #interval
    v1 = bw(p)
    v2 = bw(p + z)
    return (v2 - v1)/z

def dmiuwp (p):
    z = 0.001 #interval
    v1 = miuw(p)
    v2 = miuw(p + z)
    return (v2 - v1)/z

def drho_op (p):
    z = 0.001 #interval
    v1 = rho_oil (p)
    v2 = rho_oil (p + z)
    return (v2 - v1)/z

def drho_gp (p):
    z = 0.001 #interval
    v1 = rho_gas (p)
    v2 = rho_gas (p + z)
    return (v2 - v1)/z

def drho_wp (p):
    z = 0.001 #interval
    v1 = rho_water (p)
    v2 = rho_water (p + z)
    return (v2 - v1)/z

def drsp(p):
    z = 0.001 #interval
    v1 = rs (p)
    v2 = rs (p + z)
    return (v2 - v1)/z


# In[268]:


#fungsi subroutine rock property
def phi (p):
    phii = phi0 * np.exp (crock * (p-p_ref)) 
    return phii

def kroil (SW):
    return interpolate(Sw, Kro, SW)

def krwater (SW):
    return interpolate(Sw, Krw, SW)

def pcow (SW):
    return interpolate(Sw, Pcow, SW)

#fungsi subroutine turunan rock property
def dphip (p):
    z = 0.001 #interval
    v1 = phi (p)
    v2 = phi (p + z)
    return (v2 - v1)/z

def dkrosw (SW):
    z = 0.001 #interval
    v1 = kroil (SW)
    v2 = kroil (SW + z)
    return (v2 - v1)/z

def dkrwsw (SW):
    z = 0.001 #interval
    v1 = krwater (SW)
    v2 = krwater (SW + z)
    return (v2 - v1)/z

def dpcowsw (SW):
    z = 0.001 #interval
    v1 = pcow (SW)
    v2 = pcow (SW + z)
    return (v2 - v1)/z


# In[269]:


#fungsi subroutine perhitungan OOIP, OGIP, OWIP
def initialcond():
    global ooip, ogip, owip, P3d, S3d, dltx, dlty, dltz, tpx, tpy, tpz, V_bulk
    #interval grid
    dltx = tx / nx
    dlty = ty / ny
    dltz = tz / nz
    V_bulk = dltx * dlty * dltz
    
    tpx = 6.3283 * (10**(-3)) * kx * dlty * dltz / dltx
    tpy = 6.3283 * (10**(-3)) * ky * dltx * dltz / dlty
    tpz = 6.3283 * (10**(-3)) * kz * dltx * dlty / dltz
    
    #penentuan tekanan masing-masing kedalaman
    Pz = []
    Pz.append(Pi)
    i = 0
    while i < nz-1:
        err = 100
        Pa = Pz[i]
        Pb = Pa
        while err > 0.000001:
            Pmid = 0.5 * (Pa + Pb)
            dp = 0.001
            fPb = Pa + rho_oil (Pmid) * dltz - Pb
            Pmid2 = 0.5 * (Pa+Pb+dp)
            fPb2 = Pa + rho_oil (Pmid2)*dltz - (Pb+dp)
            dfPb = (fPb-fPb2)/(-dp)
            Pbnew = Pb - fPb / dfPb
            err = abs(Pbnew - Pb)
            Pb = Pbnew
        Pz.append (Pb)
        i = i + 1
    print(Pz)
    
    #penentuan OOIP, OGIP, dan OWIP
    P3d = np.zeros((nx, ny, nz), dtype=float) #array pressure 3 dimensi
    S3d = np.zeros((nx, ny, nz), dtype=float) #array saturasi 3 dimensi
    
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                P3d[i][j][k] = Pz[k]
                S3d[i][j][k] = Swi
    sumoil = 0
    sumgas = 0
    sumwater = 0
    
    for i in range (0, nx):
        for j in range (0, ny):
            for k in range (0, nz):
                sumoil += V_bulk * phi (P3d[i][j][k]) * (1-S3d[i][j][k]) / bo (P3d[i][j][k]) #dalam cuft
                sumgas += V_bulk * rs (P3d[i][j][k]) * phi (P3d[i][j][k]) * (1-S3d[i][j][k]) / bo (P3d[i][j][k]) #dalam cuft
                sumwater += V_bulk * phi (P3d[i][j][k]) * (S3d[i][j][k]) / bw (P3d[i][j][k]) #dalam cuft
    ooip = sumoil
    ogip = sumgas
    owip = sumwater
    #output OOIP, OGIP, dan OWIP
    print ("OOIP (MMSTB) = ", round (sumoil/(5.6146*1000000), 3))
    print ("OGIP (BSCF)  = ", round (sumgas/1000000000, 3))
    print ("OWIP (MMSTB) = ", round (sumwater/(5.6146*1000000), 3))
    return


# In[270]:


def calc_rem():
    global roip, rwip
    sumoil = 0
    sumwater = 0
    
    for i in range (0, nx):
        for j in range (0,ny):
            for k in range (0,nz):
                sumoil += V_bulk*phi(P3d[i][j][k])*(1-S3d[i][j][k]) / bo(P3d[i][j][k])
                sumwater += V_bulk*phi(P3d[i][j][k])*(1-S3d[i][j][k]) / bw(P3d[i][j][k])
    roip = sumoil
    rwip = sumwater
    return


# In[271]:


def well(t):
    global qo, qw, dQodp, dQwdp, dQods, dQwds
    
    #initialize
    
    qo = np.zeros(Nw, dtype=float)
    qw = np.zeros(Nw, dtype=float)
    dQodp = np.zeros(Nw, dtype=float)
    dQwdp = np.zeros(Nw, dtype=float)
    dQods = np.zeros(Nw, dtype=float)
    dQwds = np.zeros(Nw, dtype=float)
    
    for i in range(0, Nw):
        pwell = P3d[wlx[i]][wly[i]][wlz[i]]
        swwell = S3d[wlx[i]][wly[i]][wlz[i]]
        
        nr = 0
        while t > wr[i][nr][0]:
            nr += 1
        
        qtot = wr[i][nr][1]*5.6146
        
        if i == 0:
            wc = 1
        else:
            m = (miuo(pwell)*bo(pwell)*krwater(swwell)) / (miuw(pwell) * bw(pwell) * kroil(swwell))
            wc = m / (1+m)
        
        
        qw[i] = qtot * wc
        qo[i] = qtot * (1-wc)
        dQwds[i] = qw[i] * (1-wc) * (dkrwsw(swwell) / krwater(swwell) - dkrosw(swwell) / kroil(swwell))
        dQods[i] = -dQwds[i]
        dQwdp[i] = qw[i] * (1-wc) * (dmiuop(pwell) / miuo(pwell) + dbop(pwell) / bo(pwell) - dmiuwp(pwell) / miuw(pwell) - dbwp(pwell) / bw(pwell))
        dQodp[i] = -dQwdp[i]
    return


# In[272]:


#fungsi subroutine poten
def poten():
    global ibw, icw, idw, iew, ifw, igw, ifo, igo
    
    ibw = np.zeros((nx, ny, nz), dtype=float)
    icw = np.zeros((nx, ny, nz), dtype=float)
    idw = np.zeros((nx, ny, nz), dtype=float)
    iew = np.zeros((nx, ny, nz), dtype=float)
    ifw = np.zeros((nx, ny, nz), dtype=float)
    igw = np.zeros((nx, ny, nz), dtype=float)
    ifo = np.zeros((nx, ny, nz), dtype=float)
    igo = np.zeros((nx, ny, nz), dtype=float)
    
    for k in range (0,nz):
        for j in range (0, ny):
            for i in range (0, nx):
                ps = P3d[i][j][k] #definisi tekanan grid (diri sendiri)
                
                #ibw -> left
                if i>0:
                    pn = P3d[i-1][j][k] #definisi tekanan neighbor
                    potw = pn - ps #definisi tekanan potential
                    if potw>0:
                        ibw[i][j][k] = 1
                
                #icw -> right
                if i<(nx-1):
                    pn = P3d[i+1][j][k]
                    potw = pn - ps
                    if potw>0:
                        icw[i][j][k] = 1
                
                #idw -> back
                if j>0:
                    pn = P3d[i][j-1][k]
                    potw = pn - ps
                    if potw>0:
                        idw[i][j][k] = 1
                
                #iew -> front
                if j<(ny-1):
                    pn = P3d[i][j+1][k]
                    potw = pn - ps
                    if potw>0:
                        iew[i][j][k] = 1
                
                #ifw dan ifo -> top
                if k>0:
                    pn = P3d[i][j][k-1]
                    pm = 0.5*(pn + ps)
                    potw = pn - ps + rho_water(pm) * dltz
                    poto = pn - ps + rho_oil(pm) * dltz
                    if potw>0:
                        ifw[i][j][k] = 1
                    if poto>0:
                        ifo[i][j][k] = 1
                
                #igw dan igo -> top
                if k<(nz-1):
                    pn = P3d[i][j][k+1]
                    pm = 0.5*(pn + ps)
                    potw = pn - ps - rho_water(pm) * dltz
                    poto = pn - ps - rho_oil(pm) * dltz
                    if potw>0:
                        igw[i][j][k] = 1
                    if poto>0:
                        igo[i][j][k] = 1
    return


# In[273]:


def deriv (sws, swn, ps, pn, dz, ijkw, ijko, pgeo):
    global dT, fTw, fTo
    
    pmid = 0.5 * (pn + ps)
    if ijkw == 1:
        krwn = krwater(swn)
        dkrwn = dkrwsw(swn)
    else:
        krwn = krwater(sws)
        dkrwn = 0
    
    if ijko == 1:
        kron = kroil(swn)
        dkron = dkrosw(swn)
    else:
        kron = kroil(sws)
        dkron = 0
    
    if (dz == 0):
        dno = 0
        dnw = 0
        ddno = 0
        ddnw = 0
    else:
        dno = rho_oil(pmid)*dz
        dnw = rho_water(pmid)*dz
        ddno = drho_op(pmid)/2*dz
        ddnw = drho_wp(pmid)/2*dz
    
    Tw = krwn / (miuw(pmid)*bw(pmid))*pgeo
    To = kron / (miuo(pmid)*bo(pmid))*pgeo
    
    dT = np.zeros(4, dtype=float)
    dT[0] = dkron*pgeo*(pn-ps-dno)/(miuo(pmid)*bo(pmid))
    dT[1] = -To/2*(dmiuop(pmid)/miuo(pmid)+dbop(pmid)/bo(pmid))*(pn-ps-dno)+To*(1-ddno)
    dT[2] = dkrwn*pgeo/(miuw(pmid)*bw(pmid))*(pn-ps-dnw)
    dT[3] = -Tw/2*(dmiuwp(pmid)/miuw(pmid)+dbwp(pmid)/bw(pmid))*(pn-ps-dnw)+Tw*(1-ddnw)
    
    fTo = To*(pn-ps-dno)
    fTw = Tw*(pn-ps-dnw)
    return


# In[274]:


def jacob():
    global Ja, Jb, Jc, Jd, Je, Jf, Jg, Fw, Fo
    Ja = np.zeros((nx, ny, nz, 4), dtype=float)
    Jb = np.zeros((nx+1, ny, nz, 4), dtype=float)
    Jc = np.zeros((nx+1, ny, nz, 4), dtype=float)
    Jd = np.zeros((nx, ny+1, nz, 4), dtype=float)
    Je = np.zeros((nx, ny+1, nz, 4), dtype=float)
    Jf = np.zeros((nx, ny, nz+1, 4), dtype=float)
    Jg = np.zeros((nx, ny, nz+1, 4), dtype=float)

    fluxo = np.zeros((nx, ny, nz), dtype=float)
    fluxw = np.zeros((nx, ny, nz), dtype=float)
    
    Fo = np.zeros((nx, ny, nz), dtype=float)
    Fw = np.zeros((nx, ny, nz), dtype=float)

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                ps = P3d[i][j][k]
                sws = S3d[i][j][k]

                #untuk nb
                if i == 0:
                    bfo = 0
                    bfw = 0
                else:
                    pn = P3d[i-1][j][k]
                    swn = S3d[i-1][j][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, ibw[i][j][k], ibw[i][j][k], tpx)
                    bfo = fTo
                    bfw = fTw
                    for l in range(0, 4):
                        Jb[i][j][k][l] = dT[l]

                #untuk nc
                if i == nx-1:
                    cfo = 0
                    cfw = 0
                else:
                    pn = P3d[i+1][j][k]
                    swn = S3d[i+1][j][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, icw[i][j][k], icw[i][j][k], tpx)
                    cfo = fTo
                    cfw = fTw
                    for l in range(0, 4):
                        Jc[i][j][k][l] = dT[l]

                #untuk nd
                if j == 0:
                    dfo = 0
                    dfw = 0
                else:
                    pn = P3d[i][j-1][k]
                    swn = S3d[i][j-1][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, idw[i][j][k], idw[i][j][k], tpy)
                    dfo = fTo
                    dfw = fTw
                    for l in range(0, 4):
                        Jd[i][j][k][l] = dT[l]
                
                #untuk ne
                if j == ny-1:
                    efo = 0
                    efw = 0
                else:
                    pn = P3d[i][j+1][k]
                    swn = S3d[i][j+1][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, iew[i][j][k], iew[i][j][k], tpy)
                    efo = fTo
                    efw = fTw
                    for l in range(0, 4):
                        Je[i][j][k][l] = dT[l]
                
                #untuk nf
                if k == 0:
                    ffo = 0
                    ffw = 0
                else:
                    pn = P3d[i][j][k-1]
                    swn = S3d[i][j][k-1]
                    d = -dltz
                    deriv(sws, swn, ps, pn, d, ifw[i][j][k], ifo[i][j][k], tpz)
                    ffo = fTo
                    ffw = fTw
                    for l in range(0, 4):
                        Jf[i][j][k][l] = dT[l]
                
                # ng
                if k == nz-1:
                    gfo = 0
                    gfw = 0
                else:
                    pn = P3d[i][j][k+1]
                    swn = S3d[i][j][k+1]
                    d = dltz
                    deriv(sws, swn, ps, pn, d, igw[i][j][k], igo[i][j][k], tpz)
                    gfo = fTo
                    gfw = fTw
                    for l in range(0, 4):
                        Jg[i][j][k][l] = dT[l]

                fluxw[i][j][k] = bfw+cfw+dfw+efw+ffw+gfw
                fluxo[i][j][k] = bfo+cfo+dfo+efo+ffo+gfo

    acc = np.zeros(4, dtype=float)

    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                ps = P3d[i][j][k]
                pn = P[i][j][k]     # p previous
                sws = S3d[i][j][k]
                swn = S[i][j][k]   # sw previous

                accw = -V_bulk/dt * (phi(ps)*sws/bw(ps)-phi(pn)*swn/bw(pn))
                acco = -V_bulk/dt * (phi(ps)*(1-sws)/bo(ps)-phi(pn)*(1-swn)/bo(pn))

                acc[0] = V_bulk/dt * (phi(ps)/bo(ps))
                acc[1] = -V_bulk/dt* (1-sws)*((dphip(ps)*bo(ps)-phi(ps)*dbop(ps))/(bo(ps)**2))
                acc[2] = -V_bulk/dt*(phi(ps)/bw(ps))
                acc[3] = -V_bulk/dt * sws * ((dphip(ps)*bw(ps)-phi(ps)*dbwp(ps))/(bw(ps)**2))

                sso = 0
                ssw = 0
                ss = np.zeros(4, dtype=float)
                for l in range(0, Nw):
                    if (i == wlx[l] and j == wly[l] and k == wlz[l]):
                        sso = qo[l]
                        ssw = qw[l]
                        ss[0] = dQods[l]
                        ss[1] = dQodp[l]
                        ss[2] = dQwds[l]
                        ss[3] = dQwdp[l]
                
                for l in range(0, 4):
                    Ja[i][j][k][l] = -Jb[i+1][j][k][l] - Jc[i-1][j][k][l] - Jd[i][j+1][k][l]
                    Ja[i][j][k][l] = Ja[i][j][k][l] - Je[i][j-1][k][l] - Jf[i][j][k+1][l] - Jg[i][j][k-1][l]
                    Ja[i][j][k][l] = Ja[i][j][k][l] + acc[l]-ss[l]
                
                Fo[i][j][k] = fluxo[i][j][k] + acco - sso
                Fw[i][j][k] = fluxw[i][j][k] + accw - ssw
    return


# In[275]:


def jm_positioner():
    global jp

    # Absensi Jacobian
    jp = np.zeros((nx, ny, nz, 7), dtype=int)
    count = 0
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # Location A relative to i,j,k
                jp[i][j][k][0] = count
                count +=1
    # Neighbour Coordinator
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # Location F relative of i,j,k
                # Up
                if(k!=0):
                    jp[i][j][k][1] = jp[i][j][k-1][0]
                else:
                    jp[i][j][k][1] = -1
                # Location D relative of i,j,k
                # Back
                if(j!=0):
                    jp[i][j][k][2] = jp[i][j-1][k][0]
                else:
                    jp[i][j][k][2] = -1
                # Location B relative of i,j,k
                # Left
                if(i!=0):
                    jp[i][j][k][3] = jp[i-1][j][k][0]
                else:
                    jp[i][j][k][3] = -1
                # Location C relative of i,j,k
                # Right
                if(i!=nx-1):
                    jp[i][j][k][4] = jp[i+1][j][k][0]
                else:
                    jp[i][j][k][4] = -1
                # Location E relative of i,j,k
                # Front
                if(j!=ny-1):
                    jp[i][j][k][5] = jp[i][j+1][k][0]
                else:
                    jp[i][j][k][5] = -1
                # Location G relative of i,j,k
                # Down
                if(k!=nz-1):
                    jp[i][j][k][6] = jp[i][j][k+1][0]
                else:
                    jp[i][j][k][6] = -1
    return


# In[276]:


def jm_constructor(nx, ny, nz, Ja, Jb, Jc, Jd, Je, Jf, Jg, Fo, Fw):
    global jm, jmm
    jm = np.zeros((nx*ny*nz*2, 2*nx*ny*nz), dtype=float)
    n = 0
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                # 2-Rows per grid
                for h in range(0, 2):   # h={0,1}
                    # 7 Derivate Members
                    for m in range(0, 7):
                        # if(jp[i][j][k][m]!=-1):
                        for mm in range(0, 2):
                            if(m==0 and jp[i][j][k][m]!=-1):   # A
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 1111
                                jm[n][jp[i][j][k][m] * 2 + mm] = Ja[i][j][k][h * 2+mm]
                            elif(m==1 and jp[i][j][k][m]!=-1): # F
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 6666
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jf[i][j][k][h * 2+mm]
                            elif(m==2 and jp[i][j][k][m]!=-1): # D
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 4444
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jd[i][j][k][h * 2+mm]
                            elif(m==3 and jp[i][j][k][m]!=-1): # B
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 2222
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jb[i][j][k][h * 2+mm]
                            elif(m==4 and jp[i][j][k][m]!=-1): # C
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 3333
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jc[i][j][k][h * 2+mm]
                            elif(m==5 and jp[i][j][k][m]!=-1): # E
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 5555
                                jm[n][jp[i][j][k][m] * 2 + mm] = Je[i][j][k][h * 2+mm]
                            elif(m==6 and jp[i][j][k][m]!=-1):  # G
                                # jm[n][jp[i][j][k][m] * 2 + mm] = 7777
                                jm[n][jp[i][j][k][m] * 2 + mm] = Jg[i][j][k][h * 2+mm]
                    n+=1
    nrow = 0
    jmm = np.zeros(2*nx*ny*nz, dtype=float)
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                jmm[nrow] = -Fo[i][j][k]
                nrow+=1
                jmm[nrow] = -Fw[i][j][k]
                nrow+=1
    return


# In[277]:


# Main Program

# Collecting Arrays
aTIME = []
aDT = []
aWATINJ = []
aOILPROD = []
aWATPROD = []
aWC = []
aWOR = []
aCUMINJ = []
aCUMOPROD = []
aCUMWPROD = []
aPWBINJ = []
aPWBPROD = []
aMB_ERR_OIL = []
aMB_ERR_WAT = []

print("Subprogram:Readdata/running")
readdata()
print("Subprogram:Readdata/success")
print("")

print("Subprogram:Initial_Cond/running")
initialcond()
print("Subprogram:Initial_Cond/success")
print("")

print("Subprogram:jm_positioner/running")
jm_positioner()
print("Subprogram:jm_positioner/success")
print("")

E_s = float(0.001)  # dSw -> dianggap konvergen
E_p = float(0.1)    # dP -> dianggap konvergen
E_fo = float(1)
E_fw = float(5)

dSLIM = 0.02
dPLIM = 50

t = 0
dt = 3
tmax = 2000
cum_oilprod = 0
cum_watprod = 0
cum_watinj = 0
cum_oilinj = 0

while t<tmax:
    P = np.zeros((nx, ny, nz), dtype=float)
    S = np.zeros((nx, ny, nz), dtype=float)
    
    t = t + dt
    for k in range(0, nz):
        for j in range(0, ny):
            for i in range(0, nx):
                P[i][j][k] = P3d[i][j][k]
                S[i][j][k] = S3d[i][j][k]
    
    print("Subprogram:Poten/running")
    poten()
    print("Subprogram:Poten/success")
    print("")

    c = 0
    niter = 0
    itermax = 100
    while c == 0:
        niter += 1
        print("Subprogram:Well/running")
        well(t)
        print("Subprogram:Well/success")
        print("")

        print("Subprogram:Jacob/running")
        jacob()
        print("Subprogram:Jacob/success")
        print("")

        print("Subprogram:jm_constructor/running")
        jm_constructor(nx, ny, nz, Ja, Jb, Jc, Jd, Je, Jf, Jg, Fo, Fw)
        print("Subprogram:jm_constructor/success")
        print("")

        print("Subprogram:gauss/running")
        sol = solve(jm, jmm)  # gauss
        print("time: ", t)
        print("iter: ", niter)
        print("Subprogram:gauss/success")
        print("")

        # Update Values
        # Separate Solution to Sw & P
        x_dsw = np.zeros((nx, ny, nz), dtype=float)
        x_dp = np.zeros((nx, ny, nz), dtype=float)
        dr = 0
        for k in range(0, nz):
            for j in range(0, ny):
                for i in range(0, nx):
                    x_dsw[i][j][k] = sol[dr]
                    dr += 1
                    x_dp[i][j][k] = sol[dr]
                    dr += 1
        for k in range(0, nz):
            for j in range(0, ny):
                for i in range(0, nx):
                    S3d[i][j][k] = S3d[i][j][k]+x_dsw[i][j][k]
                    P3d[i][j][k] = P3d[i][j][k]+x_dp[i][j][k]
        x_dsw_max = np.amax(abs(x_dsw))
        x_dp_max = np.amax(abs(x_dp))
        fo_max = np.amax(abs(Fo))
        fw_max = np.amax(abs(Fw))

        if(fo_max < E_fo and fw_max < E_fw and x_dp_max < E_p and x_dsw_max < E_s):
            c = 1
        else:
            if(niter > itermax):
                dt = dt*0.5
                t = t-dt

    
    print("Subprogram:Calc_Rem/running")
    calc_rem()
    print("Subprogram:Calc_Rem/success")
    print("")

    for i in range(0, Nw):
        if qw[i] > 0:
            Qw = qw[i]
            Qo = qo[i]
            cum_watprod += Qw*dt
            cum_oilprod += Qo*dt
        if qw[i] < 0:
            Qi = abs(qw[i])
            cum_watinj += abs(qw[i])*dt
            cum_oilinj += abs(qo[i])*dt

    mbew = (owip - rwip - cum_watprod + cum_watinj)/owip
    mbeo = (ooip - roip - cum_oilprod + cum_oilinj)/ooip
    
    watcut = Qw/(Qo+Qw)
    wor = Qw/Qo

    aTIME.append(t)
    aDT.append(dt)
    aWATINJ.append(Qi)
    aOILPROD.append(Qo)
    aWATPROD.append(Qw)
    aWC.append(watcut)
    aWOR.append(wor)
    aCUMINJ.append(cum_watinj/1000)
    aCUMOPROD.append(cum_oilprod/1000)
    aCUMWPROD.append(cum_watprod/1000)
    aPWBINJ.append(P3d[0][0][4])
    aPWBPROD.append(P3d[4][4][4])
    aMB_ERR_OIL.append(mbeo)
    aMB_ERR_WAT.append(mbew)

    dPMAX = np.amax(abs(P3d-P))
    dSMAX = np.amax(abs(S3d-S))

    dtold = dt
    dT_new_p = dPLIM/dPMAX
    dT_new_s = dSLIM/dSMAX
    dt = dt*min([dT_new_s, dT_new_p])
    if(dt/dtold > 2):
        dt = dtold*2
    if(dt>30):
        dt = 30
    if(t<tmax and t+dt > tmax):
        dt = tmax - t

with open("HasilSimulasi.txt", "w+") as ww:
    ww.write("TIME DT WATINJ OILPROD WATPROD WC WOR CUMINJ CUMOPROD CUMWPROD PWBINJ PWBPROD MBEO MBEW")
    ww.write("\n")
    ww.write("Days Days SCF/D SCF/D SCF/D % SCF/SCF MSCF MSCF MSCF psia psia dec. dec.")
    ww.write("\n")
    for x in range(0, len(aTIME)):
        ww.write(str(aTIME[x])+ " ")
        ww.write(str(aDT[x])+ " ")
        ww.write(str(aWATINJ[x])+ " ")
        ww.write(str(aOILPROD[x])+ " ")
        ww.write(str(aWATPROD[x])+ " ")
        ww.write(str(aWC[x])+ " ")
        ww.write(str(aWOR[x])+ " ")
        ww.write(str(aCUMINJ[x])+ " ")
        ww.write(str(aCUMOPROD[x])+ " ")
        ww.write(str(aCUMWPROD[x])+ " ")
        ww.write(str(aPWBINJ[x])+ " ")
        ww.write(str(aPWBPROD[x])+ " ")
        ww.write(str(aMB_ERR_OIL[x])+ " ")
        ww.write(str(aMB_ERR_WAT[x])+ " ")
        ww.write("\n")

print("Reservoir Simulation Completed")


# In[279]:





# In[ ]:




