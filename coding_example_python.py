# The purpose of this script is to compare the dissociation/recombination rate coefficients between RVT and VT models

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt 
import sympy as sym
from scipy.optimize import curve_fit

def func(x,a,b,c):
        return a*x**2+b*x+c

# some useful constants
kb = 1.38e-23
cm2eV = 1.2398e-4
eV2j = 1.602e-19
cm2j = cm2eV *eV2j
kb_cm = kb/cm2j

mA = 2.65509e-26
mBC  = 2*mA

mu = mA*mBC/(mA+mBC)

V=1
h= 6.62607e-34
Edcm = 42045.0
Ed = 42045.0*cm2j

# crs: cross-section database
# Ei : rovibrational energy grids
# Et : translational energy grids
# gi : degeneracy grids
# T : temperature

# compute the overall rate coefficient at temperature T based on the cross section database
def crs_to_k(crs,Ei,Et,gi,T):
    ki=[]
    Et = np.asarray(Et)
    crs = np.asarray(crs)
    T = np.asarray(T)
    Ei = np.asarray(Ei)
    gi = np.asarray(gi)

    for i in range(Ei.size):
        # compute ki from 0 to Emax
        ki.append(compute_ki(crs[i],Et,T))
        # compute ki from Emin to Emax
        #ki.append(compute_ki_from_emin(crs[i],Et,T,Edcm - Ei[i]))


    ki = np.asarray(ki)
    k = weight_sum_ki(ki,Ei,gi,T)
    return k

# compute the state-specific rate coefficient at temperature T based on the cross section database
def compute_ki(crs,Et,T):
    C = 8*np.pi/np.sqrt(mu)*(2*np.pi*kb*T)**-1.5*1e-20
   
    Et = np.squeeze(Et)*eV2j
    

    y = np.multiply(Et,crs)
    y = np.multiply(y,(np.exp(-Et/kb/T)))


    f = integrate.simps(y, Et)

    return C*f

# compute the state-specific rate coefficient at temperature T based on the cross section database
# The integration of cross-sections starts from Emin, instead of zero energy

def compute_ki_from_emin(crs,Et,T,Emin):
    C = 8*np.pi/np.sqrt(mu)*(2*np.pi*kb*T)**-1.5*1e-20
   
    Et = np.squeeze(Et)*eV2j
    

    y = np.multiply(Et,crs)
    y = np.multiply(y,(np.exp(-Et/kb/T)))

    mask = Et<Emin*cm2j
    y[mask]=0

    f = integrate.simps(y, Et)

    return C*f

# weighted summing the state-specific rate coefficients based on the Boltzmann distribution at temperature T    
def weight_sum_ki(ki,Ei,gi,T):
    
    w = np.multiply(gi,np.exp(-Ei*cm2j/kb/T))
    wsum = np.sum(w)

    ans = np.dot(w,ki)/wsum

    return ans

# compute the overall equilibrium constant at temperature T based on rovibrational coupled database
def compute_keq(gi,Ei,T):
    gi = np.asarray(gi)
    Ei = np.asarray(Ei)*cm2j
    Qtr_O = V*(2*np.pi*mA*kb*T/h**2)**1.5
    Qel_O = 5+3*np.exp(-228/T)+np.exp(-326./T)+5*np.exp(-22800./T)+np.exp(-48600./T)
    Qo=Qtr_O*Qel_O

    Qtr_O2 = V*(2*np.pi*mBC*kb*T/h**2)**1.5
    Qrovib_O2= np.dot(gi,np.exp(-Ei/kb/T))

    Qel_O2=3+2*np.exp(-11390/T)+np.exp(-18990/T)
    Qo2=Qtr_O2*Qrovib_O2*Qel_O2

    return Qo**2/Qo2*np.exp(-(Ed/kb/T)) # diss-rec

# compute the overall equilibrium constant at temperature T based on vibrational database and rigid motor model for rotational mode
def compute_keq_decoup(gi,Ei,T):
    gi = np.asarray(gi)
    Ei = np.asarray(Ei)*cm2j
    Qtr_O = V*(2*np.pi*mA*kb*T/h**2)**1.5
    Qel_O = 5+3*np.exp(-228/T)+np.exp(-326./T)+5*np.exp(-22800./T)+np.exp(-48600./T)
    Qo=Qtr_O*Qel_O

    Qtr_O2 = V*(2*np.pi*mBC*kb*T/h**2)**1.5
    Qvib_O2= np.dot(gi,np.exp(-Ei/kb/T))
    Qrot_O2 = T /2/2.37

    Qel_O2=3+2*np.exp(-11390/T)+np.exp(-18990/T)
    Qo2=Qtr_O2*Qrot_O2*Qvib_O2*Qel_O2

    return Qo**2/Qo2*np.exp(-(Ed/kb/T)) # diss-rec    


if __name__ == '__main__':
    
    nVT=47
    nRVT= 3110
    nEtr = 32
    TList = range(1000,21000,1000)
    TList = np.asarray(TList)

    crsD_rvt = []
    crsD_vt = [] 
    crsD_vt_esp=[]
    ei_rvt = []
    ei_vt =[]

    gi_rvt =[]
    gi_vt=[]

    #####################################
    # IMPORTANT!
    # set the directory of databases used for rvt and vt models
    #####################################

    inputDir_rvt = "/hdd/tjpan/organize_crs_database_for_DSMC/output/O2O_Andrienko_RVT_crs_merge_w_fdiss_variable_fms"
    inputDir_vt = "/hdd/tjpan/organize_crs_database_for_DSMC/output/O2O_Andrienko_VT_crs_merge_w_fdiss_variable_fms"
    inputDir_eps = "/hdd/tjpan/organize_crs_database_for_DSMC/output/O2O_Andrienko_RVT_crs_merge_w_fdiss_fms"


    # readin Ei
    # ei_rvt
    buf=[]
    filename =inputDir_rvt+"/ErvList.dat" 
    f = open(filename,"r")
    for line in f:
        buf.append(([float(x) for x in line.split()]))
    f.close() 
    
    ei_rvt = [row[1] for row in buf]
    gi_rvt = [row[0] for row in buf]

    #ei_vt
    
    for T in TList :
        buf=[]
        
        filename =inputDir_vt+"/Trot"+str(int(T/1000)-1)+"/ErvList.dat" 
    
        f = open(filename,"r")
        for line in f:
            buf.append(([float(x) for x in line.split()]))
        f.close() 
        
        ####################
        # two options to define the bin energy for vt database
        # COMMENT ONE OF THEM!
       
        #use rovibrational energy as bin energy for vt model
        ei = [row[1] for row in buf]

        #use vibrational energy as bin energy for vt model
        #ei = [row[2] for row in buf]
        ####################

        ei_vt.append(ei)

        gi = [row[0] for row in buf]
        gi_vt.append(gi)
    
    # readin Etr 
    #
    
    filename =inputDir_rvt+"/EtrList.dat" 
    f = open(filename,"r")
    et=[]
    
    for line in f:
        et.append(([float(x) for x in line.split()]))
    f.close() 
    #
    
    
    # readin RVT files 
    inputDir = inputDir_rvt
    
    for preBin in range(nRVT):
        Block = 0
        inputFile = inputDir+"/crs"+"{0:0=4d}".format(preBin)+".dat"
        f = open(inputFile,"r")
        for line in f:
            if line[0]== '#':
                Block = Block +1
            else:
                if Block ==2 :
                    line = line.strip('\n')
                    crsD_rvt.append( list(map(float, line.split())))
        
    
    # readin VT files 

    for T in TList :
        crs = []
        inputDir = inputDir_vt+"/Trot"+str(int(T/1000)-1)
        
        for preBin in range(nVT):
            Block = 0
            inputFile = inputDir+"/crs"+"{0:0=4d}".format(preBin)+".dat"
            f = open(inputFile,"r")
            for line in f:
                if line[0]== '#':
                    Block = Block +1
                else:
                    if Block ==2 :
                        line = line.strip('\n')
                        crs.append( list(map(float, line.split())))
        crsD_vt.append(crs)

    
    #readin V files
    buf=[]
    filename ="/hdd/tjpan/tjpan_dissertation/RVT_VT_comparison/data/ErvList_v.dat" 
    f = open(filename,"r")
    for line in f:
        buf.append(([float(x) for x in line.split()]))
    f.close() 
    
    ei_decoup = [row[2] for row in buf]
    gi_decoup = [row[0] for row in buf]


    # compute kd and kr
    kd_rvt=[0 for i in range(len(TList))]
    kr_rvt=[0 for i in range(len(TList))]
    kd_vt=[0 for i in range(len(TList))]
    kr_vt=[0 for i in range(len(TList))]
    keq_rvt=[0 for i in range(len(TList))]
    keq_vt=[0 for i in range(len(TList))]
    keq_decoup=[0 for i in range(len(TList))]

    for idx, T in zip(range(len(TList)),TList):
        # using rvt energy levels to compute equilibrium constant, and convert every kd to kr based on the same equilibrium constant
        keq_rvt[idx] = compute_keq(gi_rvt,ei_rvt,T)
        keq_vt[idx] = compute_keq(gi_vt[idx],ei_vt[idx],T)
        keq_decoup[idx] = compute_keq_decoup(gi_decoup,ei_decoup,T)

        kd_vt[idx]=crs_to_k(crsD_vt[idx],ei_vt[idx],et,gi_vt[idx],T)
        kd_rvt[idx]=crs_to_k(crsD_rvt,ei_rvt,et,gi_rvt,T)

        kr_vt[idx]=kd_vt[idx]/keq_vt[idx]
        kr_rvt[idx]=kd_rvt[idx]/keq_rvt[idx]
############################################
    T_old = [x*1000 for x in range(1, 9)]
    kd_esp_old = [4.66386696201144e-40,	4.18489405743485e-27,	8.51167973459192e-23,	1.17774319463063e-20,	2.22866882004016e-19,	1.56362580006796e-18,	6.23051179927195e-18,	1.74427818645707e-17]
    kr_esp_old=[0 for i in range(len(T_old))]
    for idx, T in zip(range(len(T_old)),T_old):
        keq = compute_keq(gi_rvt,ei_rvt,T)
        kr_esp_old[idx] = kd_esp_old[idx]/keq

    kd_rvt_Chaithanya = [3.21519e-39, 1.3263e-26, 1.9806422e-22, 2.36346e-20, 4.10397e-19, 2.71967e-18, 1.0391e-17, 2.81306e-17]# compute from Chaithanya's code

############################################
print(keq_rvt)


############################################
# theoretical composition using different keq

Xo_rvt= [0 for i in range(len(TList))]
Xo2_rvt = [0 for i in range(len(TList))]
Xo_vt= [0 for i in range(len(TList))]
Xo2_vt = [0 for i in range(len(TList))]
Xo_decoup= [0 for i in range(len(TList))]
Xo2_decoup = [0 for i in range(len(TList))]

noi = 1e25
no2i = 1e25
ni = noi+no2i

x = sym.Symbol('x')

for i in range(len(TList)):
    
    ans = (-(keq_rvt[i]+4*noi)+np.sqrt(keq_rvt[i]**2+8*keq_rvt[i]*noi+16*keq_rvt[i]*no2i))/8

    if no2i -ans <0 :
        ans = (-(keq_rvt[i]+4*noi)-np.sqrt(keq_rvt[i]**2+8*keq_rvt[i]*noi+16*keq_rvt[i]*no2i))/8
    
    no = noi + 2*ans
    no2= no2i- ans
    n = no + no2
    Xo_rvt[i] = no/n
    Xo2_rvt[i] = no2/n

    ans = (-(keq_vt[i]+4*noi)+np.sqrt(keq_vt[i]**2+8*keq_vt[i]*noi+16*keq_vt[i]*no2i))/8

    if no2i -ans <0 :
        ans = (-(keq_vt[i]+4*noi)-np.sqrt(keq_vt[i]**2+8*keq_vt[i]*noi+16*keq_vt[i]*no2i))/8
   
    no = noi + 2*ans
    no2= no2i- ans
    n = no + no2
    Xo_vt[i] = no/n
    Xo2_vt[i] = no2/n


    ans = (-(keq_decoup[i]+4*noi)+np.sqrt(keq_decoup[i]**2+8*keq_decoup[i]*noi+16*keq_decoup[i]*no2i))/8

    if no2i -ans <0 :
        ans = (-(keq_decoup[i]+4*noi)-np.sqrt(keq_decoup[i]**2+8*keq_decoup[i]*noi+16*keq_decoup[i]*no2i))/8
   
    no = noi + 2*ans
    no2= no2i- ans
    n = no + no2
    Xo_decoup[i] = no/n
    Xo2_decoup[i] = no2/n

    

#############################################
plt.figure(1)
plt.plot(TList,kr_rvt,'bo',label='rvt')
plt.plot(TList,kr_vt,'ro',label='vt')
plt.plot(T_old,kr_esp_old,'k--',label='esp_old')

plt.xlabel('T[K]')
plt.ylabel('kr[m^6/s]')
plt.legend()
#plt.xscale('log')
#plt.show() 

f2 = plt.figure(2)
plt.plot(TList,kd_rvt,'bo',label='rvt')
plt.plot(TList,kd_vt,'ro',label='vt')
plt.plot(T_old,kd_rvt_Chaithanya,'g--',label='rvt_old')
plt.plot(T_old,kd_esp_old,'k--',label='esp_old')

plt.xlabel('T[K]')
plt.ylabel('kd[m^3/s]')
plt.legend()
plt.yscale('log')
#plt.show()       

f3 = plt.figure(3)
plt.plot(TList,keq_rvt,'bo',label='rvt')
plt.plot(TList,keq_vt,'ro',label='vt')
plt.plot(TList,keq_decoup,'go',label='decoup')


plt.xlabel('T[K]')
plt.ylabel('Keq')
plt.legend()
plt.yscale('log')

f4 = plt.figure(4)
plt.plot(TList,Xo_rvt,'bo',label='rvt')
plt.plot(TList,Xo_vt,'ro',label='vt')
plt.plot(TList,Xo_decoup,'go',label='decoup')

plt.xlabel('T[K]')
plt.ylabel('mole fraction of O')
plt.legend()

plt.show()       

