#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Riyaz Ghulam Anwary
#12217052

#TM4112
#Project Reservoir Simulation 8 2020/2021
#Assignment 8, Due date: 23/12/2020


# In[21]:


from scipy.linalg import lu_solve, lu_factor
#from scipy.linalg import solve


# In[2]:


#Read data Reservoir Simulation
import pandas as pd
import numpy as np

def read_data ():  
    global nx, ny, nz, tx, ty, tz, Pi, Swi, Phi0, Crock, Pref, Kx, Ky, Kz, pvt_gas, data, nwell, wlx, wly, wlz, nrate
    global Sw, Krw, Kro, Pcow, rhoo, rhog, rhow, Pg, Bg, Miug, pvt_oil, pvt_water, rock_prop, time, rateinj, rateprod, time_rate
    
    data = pd.read_csv('Data_ResSimNew.txt', header = None)
    pvt_oil = pd.read_csv('PVT-Oil.csv')
    pvt_water = pd.read_csv('PVT-Water.csv')

    #Number of Grid
    temp = data.iloc[5][0]
    temp = np.array(temp.split(), dtype = int)
    nx = temp[0]
    ny = temp[1]
    nz = temp[2]
    
    #dimension x, y, z in ft
    temp = data.iloc[11][0]
    temp = np.array(temp.split(), dtype = float)
    tx = temp[0]
    ty = temp[1]
    tz = temp[2]
    
    #Pi dan Swi
    temp = data.iloc[17][0]
    temp = np.array(temp.split(), dtype = float)
    Pi = temp[0]
    Swi = temp[1]
    
    #phi0, Crock, Pref
    temp = data.iloc[24][0]
    temp = np.array(temp.split(), dtype = float)
    Phi0 = temp[0]
    Crock = temp[1]
    Pref = temp [2]
    
    #Kx, Ky, Kz
    temp = data.iloc[30][0]
    temp = np.array(temp.split(), dtype = float)
    Kx = temp[0]
    Ky = temp[1]
    Kz = temp [2]
    
    #Sw, Krw, Kro, Pcow
    nrock = int(data.iloc[35][0])
    Sw = np.zeros(nrock, dtype=float)
    Krw = np.zeros(nrock, dtype=float)
    Kro = np.zeros(nrock, dtype=float)
    Pcow = np.zeros(nrock, dtype=float)
    for i in range(42, nrock + 42):
        temp = data.iloc[i][0]
        temp = np.array(temp.split(), dtype = float)
        i=i-42
        Sw[i] = temp[0]
        Krw[i] = temp[1]
        Kro[i] = temp [2]
        Pcow[i] = temp [3]
    #membuat dataframenya
    fluid_rock = {'Sw': Sw, 'Krw': Krw, 'Kro': Kro, 'Pcow': Pcow}
    rock_prop = pd.DataFrame(fluid_rock, columns = ['Sw', 'Krw', 'Kro', 'Pcow'])
    
        
    #Fluid density
    temp = data.iloc[77][0]
    temp = np.array(temp.split(), dtype = float)
    rhoo = temp[0]
    rhog = temp[1]
    rhow = temp [2]
    
    #PVTG
    npvtg = 11
    Pg = np.zeros(npvtg, dtype=float)
    Bg = np.zeros(npvtg, dtype=float)
    Miug = np.zeros(npvtg, dtype=float)
    for i in range(104, npvtg + 104):
        temp = data.iloc[i][0]
        temp = np.array(temp.split(), dtype = float)
        i=i-104
        Pg[i] = temp[0]
        Bg[i] = temp[1]
        Miug[i] = temp [2]
    #membuat dataframenya
    gas = {'P': Pg, 'Bg': Bg, 'Miug': Miug}
    pvt_gas = pd.DataFrame(gas, columns = ['P', 'Bg', 'Miug'])
    
    #NWell
    temp = data.iloc[128][0]
    temp = np.array(temp.split(), dtype = int)
    nwell = temp[0]
    
    #Well Location
    wlx = np.zeros(nwell, dtype=int)
    wly = np.zeros(nwell, dtype=int)
    wlz = np.zeros(nwell, dtype=int)
    for i in range(133, nwell + 133):
        temp = data.iloc[i][0]
        temp = np.array(temp.split(), dtype = int)
        i=i-133
        wlx[i] = temp[0]-1
        wly[i] = temp[1]-1
        wlz[i] = temp[2]-1
    
    #NRate
    temp = data.iloc[139][0]
    temp = np.array(temp.split(), dtype = int)
    nrate = temp[0]
    
    #Injector
    time = np.zeros(4, dtype=float)
    rateinj = np.zeros(4, dtype=float)
    for i in range(144, 4 + 144):
        temp = data.iloc[i][0]
        temp = np.array(temp.split(), dtype = float)
        i=i-144
        time[i] = temp[0]
        rateinj[i] = temp[1]*5.6146
    #Producer
    rateprod = np.zeros(4, dtype=float)
    for i in range(149, 4 + 149):
        temp = data.iloc[i][0]
        temp = np.array(temp.split(), dtype = float)
        i=i-149
        rateprod[i] = temp[1]*5.6146
    #membuat dataframe Injector Producer
    kol = {'t': time, 'Qinj': rateinj, 'Qprod': rateprod}
    time_rate = pd.DataFrame(kol, columns = ['t', 'Qinj', 'Qprod'])
    
    return

read_data ()


# In[3]:


def look_up(tabx, taby, x, prop):
    i = 0
    if (x<tabx[prop][i]):
        y = tabx[taby][i]
    elif (x>tabx[prop].iloc[-1]):
        y = tabx[taby].iloc[-1]
    else:
        x1 = tabx[prop][i]
        x2 = tabx[prop][i+1]
        y1 = tabx[taby][i]
        y2 = tabx[taby][i+1]
        while(x>tabx[prop][i+1]):
            i += 1
            x1 = tabx[prop][i]
            x2 = tabx[prop][i + 1]
            y1 = tabx[taby][i]
            y2 = tabx[taby][i + 1]
        y = y1 + ((y2 - y1) / (x2 - x1)) * (x - x1)
    return y


# In[4]:


#Membuat fungsi fluid property
def Bo (P): #Fungsi Formation Volume Factor
    return look_up(pvt_oil, 'Bo', P, 'P')
def Miuo (P):
    return look_up(pvt_oil, 'Miuo', P, 'P')
def Rs (P):
    Rs = look_up(pvt_oil, 'Rs', P, 'P') * 1000 / 5.6145
    return Rs
def Bw (P):
    X = pvt_water['Cw'][0] * (P - pvt_water['Pref'][0])
    Bw = pvt_water['Bw'][0] / (1 + X + ((X**2) / 2))
    return Bw
def Miuw (P):
    Y = -1 * (pvt_water['Cw'][0]) * (P - pvt_water['Pref'][0])
    Miuw = pvt_water['Miuw'][0] / (1 + Y + ((Y**2) / 2))
    return Miuw
def Rhoo (P):
    rhoo_new = ((rhoo + Rs(P) * rhog) / Bo(P)) / 144
    return rhoo_new
def Bg (P):
    return look_up(pvt_gas, 'Bg', P, 'P') * 5.615 / 1000
def Rhog (P):
    rhog_new = (rhog / Bg (P)) / 144
    return rhog_new
def Rhow (P):
    rhow_new = (rhow / Bw (P)) / 144
    return rhow_new
def Phi (P):
    phi_new = Phi0 * np.exp (Crock * (P-Pref)) 
    return phi_new
#Membuat fungsi fluid-rock property
def Kro (Sw):
    Kro = look_up(rock_prop, 'Kro', Sw, 'Sw')
    return Kro
def Krw (Sw):
    Krw = look_up(rock_prop, 'Krw', Sw, 'Sw')
    return Krw


# In[5]:


#Membuat fungsi differential fluid property
def dBop (P):
    inc = 0.0001
    dBop = (Bo (P + inc) - Bo (P)) / inc
    return dBop
def dMiuop (P):
    inc = 0.0001
    dMiuop = (Miuo (P + inc) - Miuo (P)) / inc
    return dMiuop
def dRsp (P):
    inc = 0.0001
    dRsp = (Rs (P + inc) - Rs (P)) / inc
    return dRsp
def dBwp (P):
    inc = 0.0001
    dBwp = (Bw (P + inc) - Bw (P)) / inc
    return dBwp
def dMiuwp (P):
    inc = 0.0001
    dMiuwp = (Miuw (P + inc) - Miuw (P)) / inc
    return dMiuwp
def dRhoop (P):
    inc = 0.00001
    dRhoop = (Rhoo (P + inc) - Rhoo (P)) / inc
    return dRhoop
def dRhogp (P):
    inc = 0.00001
    dRhogp = (Rhog (P + inc) - Rhog (P)) / inc
    return dRhogp
def dRhowp (P):
    inc = 0.00001
    dRhowp = (Rhow (P + inc) - Rhow (P)) / inc
    return dRhowp
def dPhip (P):
    inc = 0.00001
    dPhip = (Phi (P + inc) - Phi (P)) / inc
    return dPhip
#Membuat fungsi differential fluid-rock property
def dKrosw (Sw):
    inc = 0.00001
    dKrosw = (Kro (Sw + inc) - Kro (Sw)) / inc
    return dKrosw
def dKrwsw (Sw):
    inc = 0.00001
    dKrwsw = (Krw (Sw + inc) - Krw (Sw)) / inc
    return dKrwsw


# In[6]:


#Fungsi subroutine perhitungan OOIP, OGIP, dan OWIP
def initial():
    global sumoil, sumgas, sumwater, Pg3d, Sw3d, dltz, dlty, dltx, tpx, tpy, tpz, Vbulk
    #Besar setiap grid dan Volume setiap gridnya
    dltx = tx / nx
    dlty = ty / ny
    dltz = tz / nz
    Vbulk = dltx * dlty * dltz
    
    tpx = 6.3283 * (10**(-3)) * Kx * dlty * dltz / dltx
    tpy = 6.3283 * (10**(-3)) * Ky * dltx * dltz / dlty
    tpz = 6.3283 * (10**(-3)) * Kz * dltx * dlty / dltz
    
    #Tekanan masing-masing kedalaman
    Pz = []
    Pz.append(Pi)
    i = 0
    while i < nz-1:
        err = 100
        Pa = Pz[i]
        Pb = Pa
        while err > 0.000001: #Iterasi untuk memperoleh nilai Pb yang benar
            Pavg = 0.5 * (Pa + Pb)
            fPb = Pa + Rhoo (Pavg) * dltz - Pb
            dfPb = 0.5 * dRhoop (Pavg) * dltz - 1
            Pbnew = Pb - fPb / dfPb
            err = abs(Pbnew - Pb)
            Pb = Pbnew
        Pz.append (Pb)
        i = i + 1
    
    #Membuat array 3 dimensi tekanan dan saturasi air setiap grid
    Pg3d = np.zeros((nx, ny, nz), dtype=float)
    Sw3d = np.zeros((nx, ny, nz), dtype=float)

    #Assign nilai pada array 3 dimensi tekanan dan saturasi air setiap grid
    #Serta menghitung dan menjumlah volume minyak, gas, dan water setiap grid
    sumoil = 0
    sumgas = 0
    sumwater = 0
    for i in range (0, nx):
        for j in range (0, ny):
            for k in range (0, nz):
                Pg3d[i][j][k] = Pz[k]
                Sw3d[i][j][k] = Swi
                sumoil += Vbulk * Phi (Pg3d[i][j][k]) * (1-Sw3d[i][j][k]) / Bo (Pg3d[i][j][k])
                sumgas += Vbulk * Rs (Pg3d[i][j][k]) * Phi (Pg3d[i][j][k]) * (1-Sw3d[i][j][k]) / Bo (Pg3d[i][j][k])
                sumwater += Vbulk * Phi (Pg3d[i][j][k]) * (Sw3d[i][j][k]) / Bw (Pg3d[i][j][k])
    
    #Konversi ke satuan baku
    sumoil = round (sumoil/(5.6146*1000000), 3) #MMSTB
    sumgas = round (sumgas/1000000000, 3) #BSCF
    sumwater = round (sumwater/(5.6146*1000000), 3) #MMSTB
                
    #Mencetak output berupa OOIP, OGIP, dan OWIP
    #print ("OOIP (MMSTB): ", sumoil)
    #print ("OGIP (BSCF): ", sumgas)
    #print ("OWIP (MMSTB): ", sumwater)
    
initial()


# In[7]:


#Fungsi subroutine perhitungan potential
def poten ():
    global ibw, icw, idw, iew, ifw, igw , ifo, igo
    
    ibw = np.zeros((nx, ny, nz), dtype=float)
    icw = np.zeros((nx, ny, nz), dtype=float)
    idw = np.zeros((nx, ny, nz), dtype=float)
    iew = np.zeros((nx, ny, nz), dtype=float)
    ifw = np.zeros((nx, ny, nz), dtype=float)
    igw = np.zeros((nx, ny, nz), dtype=float)
    ifo = np.zeros((nx, ny, nz), dtype=float)
    igo = np.zeros((nx, ny, nz), dtype=float)
    
    for i in range (0, nx):
        for j in range (0, ny):
            for k in range (0, nz):
                ps = Pg3d[i][j][k]
                
                #ibw
                if i>0:
                    pn = Pg3d[i-1][j][k]
                    potw = pn - ps
                    if potw > 0:
                        ibw [i][j][k] = 1
                        
                #icw
                if i<(nx-1):
                    pn = Pg3d[i+1][j][k]
                    potw = pn - ps
                    if potw > 0:
                        icw [i][j][k] = 1
                        
                #idw
                if j>0:
                    pn = Pg3d[i][j-1][k]
                    potw = pn - ps
                    if potw > 0:
                        idw [i][j][k] = 1
                        
                #iew 
                if j<(ny-1):
                    pn = Pg3d[i][j+1][k]
                    potw = pn - ps
                    if potw > 0:
                        iew [i][j][k] = 1
                        
                #ifw, ifo
                if (k>0):
                    pn = Pg3d[i][j][k-1]
                    pm = 0.5 * (pn + ps)
                    potw = pn - ps + Rhow (pm) * dltz
                    poto = pn - ps + Rhoo (pm) * dltz
                    if potw > 0:
                        ifw[i][j][k] = 1
                    if poto > 0:
                        ifo[i][j][k] = 1
                        
                #igw, igo
                if k<(nz-1):
                    pn = Pg3d[i][j][k+1]
                    pm = 0.5 * (pn + ps)
                    potw = pn - ps - Rhow (pm) * dltz
                    poto = pn - ps - Rhoo (pm) * dltz
                    if potw > 0:
                        igw[i][j][k] = 1
                    if poto > 0:
                        igo[i][j][k] = 1
    return


# In[8]:


#Inisiasi dan Output jawaban dari fungsi potential
for i in range (0, nx):
    for j in range (0, ny):
        for k in range (0, nz):
            if i%2 == 0:
                Pg3d[i][j][k] += 10
            if j%2 == 0:
                Pg3d[i][j][k] += 10
            if k%2 == 0:
                Pg3d[i][j][k] += 10

poten()

n1 = 0
n2 = 0
n3 = 0
n4 = 0

for i in range (0, nx):
    for j in range (0, ny):
        for k in range (0, nz):
            n1 = n1 + ibw[i][j][k] + icw[i][j][k] + idw[i][j][k] + iew[i][j][k]
            n2 = n2 + ifw[i][j][k] + igw[i][j][k] + ifo[i][j][k] + igo[i][j][k]
            n3 = n3 + ibw[i][j][k] * igw[i][j][k] + icw[i][j][k] * ifw[i][j][k]
            n4 = n4 + idw[i][j][k] * igo[i][j][k] + iew[i][j][k] * ifo[i][j][k]
            
print ("n1 : " + str(n1))
print ("n2 : " + str(n2))
print ("n3 : " + str(n3))
print ("n4 : " + str(n4))
print ("ibw(1,2,3): " + str(ibw[0][1][2]))
print ("icw(2,3,4): " + str(icw[1][2][3]))
print ("idw(4,3,2): " + str(idw[3][2][1]))
print ("iew(3,2,1): " + str(iew[2][1][0]))
print ("ifw(5,1,5): " + str(ifw[4][0][4]))
print ("igw(4,2,4): " + str(ibw[3][1][3]))
print ("ifo(3,3,3): " + str(ifo[2][2][2]))
print ("igo(2,2,2): " + str(igo[1][1][1]))


# In[9]:


#Fungsi subroutine perhitungan derivative dari transmissibility
def deriv(sws, swn, ps, pn, dz, ijkw, ijko, pgeo):
    global dT, fTw, fTo
    pmid = 0.5*(pn+ps)
    if ijkw == 1:
        Krwn = Krw (swn)
        dKrwn = dKrwsw (swn)
    else :
        Krwn = Krw (sws)
        dKrwn = 0

    if ijko == 1:
        Kron = Kro (swn)
        dKron = dKrosw (swn)
    else:
        Kron = Kro (sws)
        dKron = 0

    if (dz==0):
        dno = 0
        dnw = 0
        ddno = 0
        ddnw = 0
    else:
        dno = Rhoo (pmid)*dz
        dnw = Rhow (pmid)*dz
        ddno = dRhoop (pmid)/2*dz
        ddnw = dRhowp (pmid)/2*dz
    
    Tw = Krwn/(Miuw (pmid)*Bw(pmid))*pgeo
    To = Kron/(Miuo (pmid)*Bo(pmid))*pgeo

    dT = np.zeros(4, dtype=float)
    dT[0] = dKron*pgeo*(pn-ps-dno)/(Miuo (pmid)* Bo (pmid))
    dT[1] = -To/2*(dMiuop(pmid)/Miuo(pmid)+dBop(pmid)/Bo(pmid))*(pn-ps-dno)+To*(1-ddno)
    dT[2] = dKrwn*pgeo/(Miuw(pmid)*Bw(pmid))*(pn-ps-dnw)
    dT[3] = -Tw/2*(dMiuwp(pmid)/Miuw(pmid)+dBwp(pmid)/Bw(pmid))*(pn-ps-dnw)+Tw*(1-ddnw)
    fTo = To*(pn-ps-dno)
    fTw = Tw*(pn-ps-dnw)
    return


# In[10]:


#Output Hasil Test Data dari fungsi derivative
poten()

i = 2
j = 1
k = 1

sws = Sw3d[i][j][k]+0.05
ps = Pg3d[i][j][k]

swn = Sw3d[i-1][j][k]+0.05
pn = Pg3d[i-1][j][k]
d = 0 

deriv(sws, swn, ps, pn, d, ibw[i][j][k], ibw[i][j][k], tpx)
print(dT[0], dT[1], dT[2], dT[3])
print(fTw, fTo)
sumbux = pd.DataFrame({'Sumbu X':['dT[0]', 'dT[1]', 'dT[2]', 'dT[3]','fTw', 'fTo'],
                     'Nilai':[dT[0], dT[1], dT[2], dT[3], fTw, fTo]})

swn = Sw3d[i][j-1][k]+0.05
pn = Pg3d[i][j-1][k]
d = 0 

deriv(sws, swn, ps, pn, d, idw[i][j][k], idw[i][j][k], tpy)
print(dT[0], dT[1], dT[2], dT[3])
print(fTw, fTo)
sumbuy = pd.DataFrame({'Sumbu Y':['dT[0]', 'dT[1]', 'dT[2]', 'dT[3]','fTw', 'fTo'],
                     'Nilai':[dT[0], dT[1], dT[2], dT[3], fTw, fTo]})

swn = Sw3d[i][j][k+1]+0.05
pn = Pg3d[i][j][k+1]
d = dltz 

deriv(sws, swn, ps, pn, d, igw[i][j][k], igw[i][j][k], tpz)
print(dT[0], dT[1], dT[2], dT[3])
print(fTw, fTo)
sumbuz = pd.DataFrame({'Sumbu Z':['dT[0]', 'dT[1]', 'dT[2]', 'dT[3]','fTw', 'fTo'],
                     'Nilai':[dT[0], dT[1], dT[2], dT[3], fTw, fTo]})


# In[11]:


#Fungsi Interpolasi rate
def Qprod (t):
    return look_up(time_rate, 'Qprod', t, 't')
def Qinj (t):
    return look_up(time_rate, 'Qinj', t, 't')


# In[12]:


def well(t):
    global qo, qw, dQodp, dQwdp, dQods, dQwds
    
    #inisialisasi
    
    qo = np.zeros(nwell, dtype=float)
    qw = np.zeros(nwell, dtype=float)
    dQodp = np.zeros(nwell, dtype=float)
    dQwdp = np.zeros(nwell, dtype=float)
    dQods = np.zeros(nwell, dtype=float)
    dQwds = np.zeros(nwell, dtype=float)
    
    for i in range(0, nwell):
        pwell = Pg3d[wlx[i]][wly[i]][wlz[i]]
        swwell = Sw3d[wlx[i]][wly[i]][wlz[i]]
        
        if i == 0:
            wc = 1
            qtot = Qinj (t)
        else:
            mob = (Miuo(pwell) * Bo(pwell) / Kro(swwell)) / (Miuw(pwell) * Bw(pwell) / Krw(swwell))
            wc = mob / (1 + mob)
            qtot = Qprod (t)
            
        qw[i] = qtot * wc
        qo[i] = qtot * (1-wc)
        dQwds[i] = qw[i] * (1-wc) * (dKrwsw(swwell) / Krw(swwell) - dKrosw(swwell) / Kro(swwell))
        dQods[i] = -dQwds[i]
        dQwdp[i] = qw[i] * (1-wc) * (dMiuop(pwell) / Miuo(pwell) + dBop(pwell) / Bo(pwell) - dMiuwp(pwell) / Miuw(pwell) - dBwp(pwell) / Bw(pwell))
        dQodp[i] = -dQwdp[i]


# In[13]:


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
                ps = Pg3d[i][j][k]
                sws = Sw3d[i][j][k]

                # nb
                if i==0:
                    bfo = 0
                    bfw = 0
                else:
                    pn = Pg3d[i-1][j][k]
                    swn = Sw3d[i-1][j][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, ibw[i][j][k], ibw[i][j][k], tpx)
                    bfo = fTo
                    bfw = fTw
                    for l in range(0, 4):
                        Jb[i][j][k][l] = dT[l]

                # nc
                if i==nx-1:
                    cfo = 0
                    cfw = 0
                else:
                    pn = Pg3d[i+1][j][k]
                    swn = Sw3d[i+1][j][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, icw[i][j][k], icw[i][j][k], tpx)
                    cfo = fTo
                    cfw = fTw
                    for l in range(0, 4):
                        Jc[i][j][k][l] = dT[l]

                # nd
                if j==0:
                    dfo = 0
                    dfw = 0
                else:
                    pn = Pg3d[i][j-1][k]
                    swn = Sw3d[i][j-1][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, idw[i][j][k], idw[i][j][k], tpy)
                    dfo = fTo
                    dfw = fTw
                    for l in range(0, 4):
                        Jd[i][j][k][l] = dT[l]
                
                # ne
                if j==ny-1:
                    efo = 0
                    efw = 0
                else:
                    pn = Pg3d[i][j+1][k]
                    swn = Sw3d[i][j+1][k]
                    d = 0
                    deriv(sws, swn, ps, pn, d, iew[i][j][k], iew[i][j][k], tpy)
                    efo = fTo
                    efw = fTw
                    for l in range(0, 4):
                        Je[i][j][k][l] = dT[l]
                
                # nf
                if k==0:
                    ffo = 0
                    ffw = 0
                else:
                    pn = Pg3d[i][j][k-1]
                    swn = Sw3d[i][j][k-1]
                    d = -dltz
                    deriv(sws, swn, ps, pn, d, ifw[i][j][k], ifo[i][j][k], tpz)
                    ffo = fTo
                    ffw = fTw
                    for l in range(0, 4):
                        Jf[i][j][k][l] = dT[l]
                
                # ng
                if k==nz-1:
                    gfo = 0
                    gfw = 0
                else:
                    pn = Pg3d[i][j][k+1]
                    swn = Sw3d[i][j][k+1]
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
                ps = Pg3d[i][j][k]
                pn = P[i][j][k]     # p previous
                sws = Sw3d[i][j][k]
                swn = S[i][j][k]   # sw previous

                accw = -Vbulk/dt*(Phi(ps)*sws/Bw(ps)-Phi(pn)*swn/Bw(pn))
                acco = -Vbulk/dt*(Phi(ps)*(1-sws)/Bo(ps)-Phi(pn)*(1-swn)/Bo(pn))

                acc[0] = Vbulk/dt*(Phi(ps)/Bo(ps))
                acc[1] = -Vbulk/dt*(1-sws)*((dPhip(ps)*Bo(ps)-Phi(ps)*dBop(ps))/(Bo(ps)*Bo(ps)))
                acc[2] = -Vbulk/dt*(Phi(ps)/Bw(ps))
                acc[3] = -Vbulk/dt*sws*((dPhip(ps)*Bw(ps)-Phi(ps)*dBwp(ps))/(Bw(ps)**2))

                sso = 0
                ssw = 0
                ss = np.zeros(4, dtype=float)
                for l in range(0, nwell):
                    if (i==wlx[l] and j==wly[l] and k==wlz[l]):
                        sso = qo[l]
                        ssw = qw[l]
                        ss[0] = dQods[l]
                        ss[1] = dQodp[l]
                        ss[2] = dQwds[l]
                        ss[3] = dQwdp[l]
                
                for l in range(0, 4):
                    Ja[i][j][k][l] = -Jb[i+1][j][k][l]-Jc[i-1][j][k][l]-Jd[i][j+1][k][l]
                    Ja[i][j][k][l] = Ja[i][j][k][l]-Je[i][j-1][k][l]-Jf[i][j][k+1][l]-Jg[i][j][k-1][l]
                    Ja[i][j][k][l] = Ja[i][j][k][l] + acc[l]-ss[l]
                
                Fo[i][j][k] = fluxo[i][j][k]+acco-sso
                Fw[i][j][k] = fluxw[i][j][k]+accw-ssw
    return


# In[14]:


initial()
S = np.zeros((nx, ny, nz), dtype=float)
P = np.zeros((nx, ny, nz), dtype=float)
for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                S[i][j][k] = Sw3d[i][j][k] 
                P[i][j][k] = Pg3d[i][j][k]
                Sw3d[i][j][k] = S[i][j][k] + 0.01
                Pg3d[i][j][k] = P[i][j][k] - 20
                
time = 800
dt = 0.1
i = 0
j = 0
k = 4

Pg3d[i][j][k] -= 20

poten()
well(time)
jacob()

for l in range (0, 4):
    print(Jb[i][j][k][l], Jc[i][j][k][l], Jd[i][j][k][l])
    print(Je[i][j][k][l], Jf[i][j][k][l], Jg[i][j][k][l])
    print(Ja[i][j][k][l])
    print("")
    
print (Fw[i][j][k], Fo[i][j][k])
test1 = pd.DataFrame({'Test 1':['Jb', 'Jc', 'Jd', 'Je','Jf', 'Jg', 'Ja', 'Fw', 'Fo'],
                     '0':[Jb[i][j][k][0], Jc[i][j][k][0], Jd[i][j][k][0], Je[i][j][k][0], Jf[i][j][k][0], Jg[i][j][k][0], Ja[i][j][k][0], Fw[i][j][k], Fo[i][j][k]],
                     '1':[Jb[i][j][k][1], Jc[i][j][k][1], Jd[i][j][k][1], Je[i][j][k][1], Jf[i][j][k][1], Jg[i][j][k][1], Ja[i][j][k][1], Fw[i][j][k], Fo[i][j][k]],
                     '2':[Jb[i][j][k][2], Jc[i][j][k][2], Jd[i][j][k][2], Je[i][j][k][2], Jf[i][j][k][2], Jg[i][j][k][2], Ja[i][j][k][2], Fw[i][j][k], Fo[i][j][k]],
                     '3':[Jb[i][j][k][3], Jc[i][j][k][3], Jd[i][j][k][3], Je[i][j][k][3], Jf[i][j][k][3], Jg[i][j][k][3], Ja[i][j][k][3], Fw[i][j][k], Fo[i][j][k]]})
print("")

time = 800
dt = 1
i = 4
j = 4
k = 0

Pg3d[i][j][k] -= 20

poten()
well(time)
jacob()

for l in range (0, 4):
    print(Jb[i][j][k][l], Jc[i][j][k][l], Jd[i][j][k][l])
    print(Je[i][j][k][l], Jf[i][j][k][l], Jg[i][j][k][l])
    print(Ja[i][j][k][l])
    print("")

print (Fw[i][j][k], Fo[i][j][k])
test2 = pd.DataFrame({'Test 2':['Jb', 'Jc', 'Jd', 'Je','Jf', 'Jg', 'Ja', 'Fw', 'Fo'],
                     '0':[Jb[i][j][k][0], Jc[i][j][k][0], Jd[i][j][k][0], Je[i][j][k][0], Jf[i][j][k][0], Jg[i][j][k][0], Ja[i][j][k][0], Fw[i][j][k], Fo[i][j][k]],
                     '1':[Jb[i][j][k][1], Jc[i][j][k][1], Jd[i][j][k][1], Je[i][j][k][1], Jf[i][j][k][1], Jg[i][j][k][1], Ja[i][j][k][1], Fw[i][j][k], Fo[i][j][k]],
                     '2':[Jb[i][j][k][2], Jc[i][j][k][2], Jd[i][j][k][2], Je[i][j][k][2], Jf[i][j][k][2], Jg[i][j][k][2], Ja[i][j][k][2], Fw[i][j][k], Fo[i][j][k]],
                     '3':[Jb[i][j][k][3], Jc[i][j][k][3], Jd[i][j][k][3], Je[i][j][k][3], Jf[i][j][k][3], Jg[i][j][k][3], Ja[i][j][k][3], Fw[i][j][k], Fo[i][j][k]]})


# In[15]:


#Output hasil dalam dataframe dan write to csv
print(test1)
print(test2)
test1.to_csv('Test1_Output.csv', index = False)
test2.to_csv('Test2_Output.csv', index = False)


# In[29]:


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


# In[31]:


def jm_constructor2(nx, ny, nz, Ja, Jb, Jc, Jd, Je, Jf, Jg, Fo, Fw):
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
    # if nnn == 1:
    #     with open("jcb.txt", "w+") as ww:
    #         for r in range(0, 2*Ngx*Ngy*Ngz):
    #             for c in range(0, 2*Ngx*Ngy*Ngz):
    #                 if(c!=2*Ngx*Ngy*Ngz-1):
    #                     ww.write(str(jm[r][c])+" ")
    #                 else:
    #                     ww.write(str(jm[r][c]))
    #                     ww.write("\n")
    return


# In[32]:


def calc_rem():
    global roip, rwip
    sumoil = 0
    sumwat = 0

    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                sumoil += Vbulk * Phi (Pg3d[i][j][k]) * (1-Sw3d[i][j][k]) / Bo (Pg3d[i][j][k])
                sumwat += Vbulk * Phi (Pg3d[i][j][k]) * (Sw3d[i][j][k]) / Bw (Pg3d[i][j][k])
    roip = round(sumoil/5.6146, 3)
    rwip = round(sumwat/5.6146, 3)
    return


# In[33]:


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
read_data()
print("Subprogram:Readdata/success")
print("")

print("Subprogram:Initial_Cond/running")
initial()
print("Subprogram:Initial_Cond/success")
print("")

print("Subprogram:jm_positioner/running")
jm_positioner()
print("Subprogram:jm_positioner/success")
print("")


E_s = float(0.001)  # dSw untuk dianggap konvergen
E_p = float(0.1)    # dP utk dianggap konvergen
E_fo = float(1)
E_fw = float(5)

dSLIM = 0.02
dPLIM = 50


# In[ ]:


print("Subprogram:gauss/running")
#sol = solve(jm, jmm)  # gauss
lu, piv = lu_factor(jm) # LU Decomposition
sol = lu_solve((lu, piv), jmm)
print("time: ", t)
print("iter: ", niter)
print("Subprogram:gauss/success")
print("")

