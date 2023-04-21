"""
Based on the heat current auto-correlation function (HCACF) obtained by the 
"ave/correlate" command from LAMMPS to calculate the Run thermal conductivity (RTC), 
the format of HCAC file is "# Index TimeDelta Ncount c_flux[1]*c_flux[1] c_flux[2]*c_flux[2] c_flux[3]*c_flux[3]"

Created on Sun Dec 18 19:54:03 2022

@author: eastshng (eastsheng@homail.com)
"""
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
import logging

def read_vol(lammpsdata):
    with open(lammpsdata,"r") as f:
        for line in f:
            if "xlo" in line:
                xlo = float(line.split()[0])
                xhi = float(line.split()[1])
            if "ylo" in line:
                ylo = float(line.split()[0])
                yhi = float(line.split()[1])
            if "zlo" in line:
                zlo = float(line.split()[0])
                zhi = float(line.split()[1])
        a = xhi-xlo
        b = yhi-ylo
        c = zhi-zlo
        V = a*b*c
        print("x:%s \ny:%s \nz:%s \nvolume:%s" %(a,b,c,V))
    return V


class EMDGreenKuboTC(object):
    """docstring for ClassName"""
    def __init__(self, hac_file,timestep,
                    temperature,vol):
        super(EMDGreenKuboTC, self).__init__()
        self.hac_file = hac_file
        self.timestep = timestep
        self.temperature = temperature
        self.vol = vol

    def read_hac(self):
        with open(self.hac_file,'r') as f:
            lines = f.readlines()
            # # correlation length p
            self.p = int(lines[3].split()[1])
            # # sample interval
            self.s = int(lines[5].split()[1])-int(lines[4].split()[1])
            # # dump interval
            self.d = int(self.p*self.s)
            Nframe = int(lines[3].split()[0])
            Mframe = int(lines[-self.p-1].split()[0])
            self.nf = int((Mframe-Nframe)/self.d)
        print("Correlation length =",self.p,"step")
        print("sample interval =",self.s,"step")
        print("dump interval =",self.d,"step")
        print("number of frame =",self.nf)
        # # skip the initial zero values
        hac_matrix = np.zeros((self.p, self.nf))
        for i in tqdm(range(self.nf)):
            # print(i)
            hac = np.loadtxt(self.hac_file,skiprows=3+1+(self.p+1)*(i+1),max_rows=self.p)[:,-3:]
            hac_xx,hac_yy,hac_zz = hac[:,0],hac[:,1],hac[:,2]
            hac_xyz = (hac_xx + hac_yy + hac_zz)
            hac_matrix[:,i] = hac_xyz ## each column denotes each independent simulations
            # time.sleep(0.0005)
        hac_ave = np.mean(hac_matrix, axis=1).reshape(-1,1)
        t = np.arange(1,self.p+1)*self.s*self.timestep*1e-3
        t = t.reshape(-1,1) # # the correlation time ps
        hac_ave = np.hstack((t,hac_ave))
        return hac_matrix, hac_ave

    def green_kubo_tc(self,hac_matrix):
        kB = 1.3806504e-23
        kcalpermol2J = 4186.0/6.02214e23
        fs2s = 1.0e-15
        A2m  = 1.0e-10
        scale = self.s*self.timestep/kB/self.temperature/self.temperature/self.vol/3.0
        convert_kappa = kcalpermol2J*kcalpermol2J/fs2s/A2m

        rtc = scale*convert_kappa*integrate.cumtrapz(hac_matrix, axis = 0, initial=0) # obtain the running thermal conductivity for all simulations
        rtc_ave = np.mean(rtc, axis=1)

        # kappa_product = np.mean(rtc_ave[int(0.5*self.p + 1):])
        # kappa_coverged = np.mean(rtc[int(0.5*self.p + 1):], axis = 0)
        # kappa_error    = np.std(kappa_coverged)/np.sqrt(self.nf)

        # print(r"The thermal conductivity is calculated as %s (W/(mK)), error is %s\n" %(kappa_product, kappa_error))

        return rtc,rtc_ave

        
if __name__ == '__main__':
    # # ------------- variables ------------- 
    dt = 4.0      # time step, unit in fs
    T  = 70      # temperature, unit in K
    ks = 1        # the number of heat flux files

    path="./"
    m_hac,m_rtc = [],[]
    # for k in range(21,ks+1):
    for k in range(1,ks+1):
        print(k,"-----",ks)
        lammpsdata = str(k)+"_nve_product.data"
        V = read_vol(path+lammpsdata)
        # # ------------- read HCACF data ------------- 
        egtc = EMDGreenKuboTC(path+"J0Jt.dat",
            timestep = dt,temperature = T, vol = V)
        hac_matrix, hac_ave = egtc.read_hac()
        t = hac_ave[:,0].reshape(-1,1)
        # # ------------- Normalized HCACF data ------------- 
        hac_matrix_normalized = hac_matrix/hac_matrix[0]
        hac_normalized = hac_ave[:,-1]#/hac_ave[:,-1][0]
        hac_normalized = hac_normalized.tolist()
        # # ------------- calculating RTC -------------
        rtc,rtc_ave = egtc.green_kubo_tc(hac_matrix)
        rtc_ave = rtc_ave.tolist()
        m_hac.append(hac_normalized)
        m_rtc.append(rtc_ave)
    m_hac = np.array(m_hac).flatten('F').reshape(-1,ks)
    m_rtc = np.array(m_rtc).flatten('F').reshape(-1,ks)
    m,n = m_rtc.shape
    hac_ks_ave = np.mean(m_hac, axis=1)
    rtc_ks_ave = np.mean(m_rtc, axis=1)

    kappa_product  = np.mean(rtc_ks_ave[int(0.5*m + 1):])
    kappa_coverged = np.mean(m_rtc[int(0.5*m + 1):], axis = 0)
    kappa_error    = np.std(kappa_coverged)/np.sqrt(n)
    logging.basicConfig(filename=path+str(ks)+'_rtc.log', filemode='w', level=logging.INFO)
    logging.info(r"The average thermal conductivity is calculated as %s (W/(mK)), error is %s" %(kappa_product, kappa_error))
    print(r"The average thermal conductivity is calculated as %s (W/(mK)), error is %s" %(kappa_product, kappa_error))

    # # ------------------- save HCACF and RTC DATA ----------------------
    m_hac = np.hstack((t,m_hac))
    m_rtc = np.hstack((t,m_rtc))
    print(m_hac.shape,m_rtc.shape)
    np.savetxt(path+"all_hac_"+str(ks)+".dat",m_hac,fmt="%f")
    np.savetxt(path+"all_rtc_"+str(ks)+".dat",m_rtc,fmt="%f")

    # # # ------------------- plot ----------------------
    plt.rc('font', family='Times New Roman', size=20)
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(bottom=0.2,left=0.15)
    ax1=fig.add_subplot(211)
    # # Normalied HCACF
    for i in range(1,m_hac.shape[1]):
        plt.plot(t,m_hac[:,i]/m_hac[:,i][0], linewidth=0.3,color="grey",alpha=0.3)
    ax1.plot(t,hac_ks_ave/hac_ks_ave[0], color="blue", linewidth=2)
    ax1.set_ylabel('HCACF',fontweight='bold',size=26)
    # ax1.set_xlabel('Time (ps)',fontweight='bold',size=26) 
    # ax1.set_xticks([0,5,10,15,20])
    # ax1.set_xlim(0,0.5)

    ax2=fig.add_subplot(212)
    # # runing thermal conductivity (RTC)
    for i in range(1,m_rtc.shape[1]):
        plt.plot(t,m_rtc[:,i], linewidth=0.3,color="grey",alpha=0.3)    
    ax2.plot(t,rtc_ks_ave,color="red", linewidth=2)
    ax2.set_ylabel(r"$\regular\kappa$ (W/mK)",fontweight='bold',size=26)
    ax2.set_xlabel('Time (ps)',fontweight='bold',size=26) 
    # ax2.set_xticks([0,5,10,15,20])
    plt.savefig(path+"hac_rtc_"+str(ks)+".png")
    plt.show()