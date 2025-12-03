"""
Green-Kubo tc based the heat flux from thermo command
"""
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import signal
from pathlib import Path
from readlog import ReadLog
from tqdm import tqdm

# convert timelapse from sec to H:M:S
def convert_time(seconds):
    seconds = seconds % (24 * 3600)
    hour = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return "%d:%02d:%02d" % (hour, minutes, seconds)

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
        # V = a*b*(c-22)
        # V = a*b*(22)
        print("x:%s \ny:%s \nz:%s \nvolume:%s" %(a,b,c,V))
    return V



def read_heatflux(logfile,hf_labelx,hf_labely,hf_labelz,vol,nlog=3):
	rl = ReadLog(logfile) 
	pd_thermo = rl.ReadThermo(nlog)
	
	hfx = np.transpose(pd_thermo[[hf_labelx]].values)[0]#*vol
	hfy = np.transpose(pd_thermo[[hf_labely]].values)[0]#*vol
	hfz = np.transpose(pd_thermo[[hf_labelz]].values)[0]#*vol

	return hfx,hfy,hfz

def autocorr(Jx,mode):
	l = len(Jx)
	ac = signal.correlate(Jx,Jx, mode=mode)[-l:]
	acx = ac/l #np.arange(l,0,-1)	
	return acx

def hfacf(Jx,Jy,Jz,mode="full"):
	acx = autocorr(Jx, mode=mode)
	acy = autocorr(Jy, mode=mode)
	acz = autocorr(Jz, mode=mode)

	return acx, acy, acz

def Normalized(vec):
	return vec/vec[0]

# def plot_full_hfacf(Jx,Jy,Jz,thermointer,timestep):
# 	plt.rc('font', family='Times New Roman', size=20)
# 	fig = plt.figure(figsize=(8,6))
# 	fig.subplots_adjust(bottom=0.2,left=0.15)
# 	ax1=fig.add_subplot(111)

# 	acx, acy, acz = hfacf(Jx,Jy,Jz,mode="full")
# 	# acx = Normalized(acx)
# 	# acy = Normalized(acy)
# 	# acz = Normalized(acz)
# 	t = np.arange(1,len(acx)+1)*thermointer*timestep # fs
# 	ax1.plot(t*1e-3,acx, "r--", linewidth=0.5, label="in x-direction")
# 	ax1.plot(t*1e-3,acy, "g--", linewidth=0.5, label="in y-direction")
# 	ax1.plot(t*1e-3,acz, "b--", linewidth=0.5, label="in z-direction")
# 	ax1.plot(t*1e-3,(acx+acy+acz)/3, "r-", linewidth=1, label="ave")

# 	# plt.xscale('log')
# 	# plt.ylim([-.1,1.2])
# 	# plt.xlim(0.001,500)
# 	plt.legend()
# 	ax1.set_ylabel("HFACF", fontweight="bold",fontsize =22)
# 	ax1.set_xlabel("Correlation Time (ps)", fontweight="bold", fontsize = 22)
# 	# plt.savefig("HFACF.png", dpi=300)
# 	plt.show()
# 	return


# def plot_hfacf(Jx,Jy,Jz,sample_interval,correlate_step,thermointer,timestep):
# 	plt.rc('font', family='Times New Roman', size=20)
# 	fig = plt.figure(figsize=(8,6))
# 	fig.subplots_adjust(bottom=0.2,left=0.15)
# 	ax1=fig.add_subplot(111)
# 	number_of_ave = int(np.floor(len(Jx)/correlate_step/sample_interval*thermointer))	
# 	acx_s,acy_s,acz_s = 0,0,0
# 	for i in tqdm(range(number_of_ave)):
# 		Jxi = Jx[int(i*correlate_step):int((i+1)*correlate_step)]
# 		Jyi = Jy[int(i*correlate_step):int((i+1)*correlate_step)]
# 		Jzi = Jz[int(i*correlate_step):int((i+1)*correlate_step)]
# 		acx, acy, acz = hfacf(Jxi,Jyi,Jzi,mode="full")
# 		acx_s += acx
# 		acy_s += acy
# 		acz_s += acz
# 	acx = acx_s/number_of_ave	
# 	acy = acy_s/number_of_ave	
# 	acz = acz_s/number_of_ave	
# 	# acx = Normalized(acx)
# 	# acy = Normalized(acy)
# 	# acz = Normalized(acz)
# 	ac_ave = (acx+acy+acz)/3
# 	t = np.arange(1,correlate_step+1)*sample_interval*timestep # fs
# 	HFACF = np.vstack((t,acx,acy,acz,ac_ave)).T
# 	np.savetxt("hfacf.dat",HFACF,fmt="%f")
# 	ax1.plot(t*1e-3,acx, "r--", linewidth=0.5, label="in x-direction")
# 	ax1.plot(t*1e-3,acy, "g--", linewidth=0.5, label="in y-direction")
# 	ax1.plot(t*1e-3,acz, "b--", linewidth=0.5, label="in z-direction")
# 	ax1.plot(t*1e-3,ac_ave, "r-", linewidth=1, label="ave")

# 	# plt.xscale('log')
# 	# plt.ylim([-.1,1.2])
# 	# plt.xlim(0.001,500)
# 	plt.legend()
# 	ax1.set_ylabel("HFACF", fontweight="bold",fontsize =22)
# 	ax1.set_xlabel("Correlation Time (ps)", fontweight="bold", fontsize = 22)
# 	# plt.savefig("HFACF.png", dpi=300)
# 	plt.show()
# 	return

def plot_hfacf_intersect(hfacf_file,Jx,Jy,Jz,sample_interval,correlate_step,thermointer,timestep,intersect_factor=5):
	plt.rc('font', family='Times New Roman', size=20)
	fig = plt.figure(figsize=(8,6))
	fig.subplots_adjust(bottom=0.2,left=0.15)
	ax1=fig.add_subplot(111)
	acx_s = np.zeros(correlate_step)
	acy_s = np.zeros(correlate_step)
	acz_s = np.zeros(correlate_step)
	intersect = int(correlate_step/intersect_factor)
	number_of_ave = int((len(Jx)-correlate_step)/intersect)
	t = np.arange(1,correlate_step+1)*sample_interval*timestep # fs
	count = 0
	all_hf=open(hfacf_file.split(".")[0]+"_all.dat","w")
	for i in tqdm(range(number_of_ave)):
		m = int(i*(intersect))
		n = m+correlate_step
		if n<=len(Jx):
			count += 1
			# print(m,n)
			Jxi = Jx[m:n]
			Jyi = Jy[m:n]
			Jzi = Jz[m:n]
			acxi, acyi, aczi = hfacf(Jxi,Jyi,Jzi,mode="full")
			
			# ax1.plot(t*1e-3,Normalized(acxi), "r--", linewidth=0.1, alpha=0.8)
			# ax1.plot(t*1e-3,Normalized(acyi), "g--", linewidth=0.1, alpha=0.8)
			# ax1.plot(t*1e-3,Normalized(aczi), "b--", linewidth=0.1, alpha=0.8)

			all_hf.write("Frame:"+str(i)+"\n")
			for j in range(len(acxi)):
				all_hf.write(str(acxi[j])+"\t"+str(acyi[j])+"\t"+str(aczi[j])+"\n")
			
			acx_s += acxi
			acy_s += acyi
			acz_s += aczi
	# print(acy_s.shape)
	acx = acx_s/count
	acy = acy_s/count
	acz = acz_s/count

	# print(count,number_of_ave)
	ac_ave = (acx+acy+acz)/3
	HFACF = np.vstack((t,acx,acy,acz,ac_ave)).T
	np.savetxt(hfacf_file,HFACF,fmt="%f")

	ax1.plot(t*1e-3,Normalized(acx), "r-", linewidth=0.8, label=r'$\regular{H_x}$')
	ax1.plot(t*1e-3,Normalized(acy), "g-", linewidth=0.8, label=r'$\regular{H_y}$')
	ax1.plot(t*1e-3,Normalized(acz), "b-", linewidth=0.8, label=r'$\regular{H_z}$')
	ax1.plot(t*1e-3,Normalized(ac_ave), "r-", linewidth=2, label=r'$\regular{H_{ave}}$')
	# plt.xscale('log')
	# plt.ylim([-.1,1.2])
	# plt.xlim(0.001,500)
	plt.legend()
	ax1.set_ylabel("HFACF", fontweight="bold",fontsize =22)
	ax1.set_xlabel("Correlation Time (ps)", fontweight="bold", fontsize = 22)
	# plt.savefig("HFACF.png", dpi=300)
	# plt.show()
	all_hf.close()
	return


def green_kubo_tc(hfacf_file,runingtc_file,sample_interval,timestep,TK,vol,ave_factor=0.5):
	"""
	kappa = (V/kB/T/T)*integrate(<J0.Jt>) = (1/kB/T/T/V)*integrate(<cJ0.cJt>)
	"""
	HFACF = np.loadtxt(hfacf_file)
	t, acx, acy, acz, ac_ave = HFACF[:,0],HFACF[:,1],HFACF[:,2],HFACF[:,3],HFACF[:,4]

	kB = 1.3806504e-23 # J/K
	kcalpermol2J = 4186.0/6.02214e23 # kcal/mol to J
	fs2s = 1.0e-15
	A2m  = 1.0e-10
	convert_to_SI = kcalpermol2J*kcalpermol2J/fs2s/A2m
	scale = convert_to_SI/kB/TK/TK/vol*sample_interval*timestep
	rtcX = integrate.cumtrapz(acx, initial=0)*scale
	rtcY = integrate.cumtrapz(acy, initial=0)*scale
	rtcZ = integrate.cumtrapz(acz, initial=0)*scale
	# rtc = integrate.cumtrapz(ac_ave, initial=0)*scale*convert_to_SI
	rtc = (rtcX+rtcY+rtcZ)/3
	m = len(rtc)
	kappa_product  = np.mean(rtc[int(ave_factor*m + 1):])
	kappa_error    = np.std(rtc[int(ave_factor*m + 1):])
	runingtc = np.vstack((t,rtcX,rtcY,rtcZ,rtc)).T
	np.savetxt(runingtc_file,runingtc,fmt="%f",
		header="kappa = "+str(kappa_product)+" W/mK, error = "+str(kappa_error))
	
	print("kappa = %s (W/mK), error = ± %s" %(kappa_product,kappa_error))

	plt.rc('font', family='Times New Roman', size=20)
	fig = plt.figure(figsize=(8,6))
	fig.subplots_adjust(bottom=0.2,left=0.15)
	ax1=fig.add_subplot(111)     
 
	ax1.plot(t*1e-3,rtcX, "r-", linewidth=0.8, label=r'$\regular{\kappa_x}$',alpha = 0.3)
	ax1.plot(t*1e-3,rtcY, "g-", linewidth=0.8, label=r'$\regular{\kappa_y}$',alpha = 0.3)
	ax1.plot(t*1e-3,rtcZ, "b-", linewidth=0.8, label=r'$\regular{\kappa_z}$',alpha = 0.3)
	ax1.plot(t*1e-3,rtc, "r-", linewidth=2.0, label=r'$\regular{\kappa_{ave}}$')
	# plt.xscale('log')
	# plt.ylim([-.1,1.2])
	# plt.xlim(0.001,500)
	plt.legend()
	ax1.set_ylabel(r"$\regular{\kappa}$ (W/mK)", fontweight="bold",fontsize =22)
	ax1.set_xlabel("Correlation Time (ps)", fontweight="bold", fontsize = 22)
	# plt.savefig("HFACF.png", dpi=300)
	plt.show()
	return


if __name__ == '__main__':
	## ------ runing time setting ------ ##
	t_start = time.perf_counter()
	## ------ path setting ------ ##
	path = "./"
	logfile = "1_green_kubo_tc_thermo_HeatFlux.lammps"
	# Path(path+"fig/").mkdir(parents=True,exist_ok=True)
	
	## ------ lammps setting ------ ##
	timestep = 4 # fs
	thermointer = 10 # thermo inter / step
	TK  = 70
	lammpsdata = "1_nve_product.data"

	## ------ correlate setting ------ ##
	sample_interval = 10 # step
	correlate_step = 200 # step

	## ------ output setting ------ ##
	hfacf_file = "hfacf.dat"
	runingtc_file = "runingtc.dat"

	## ------ green-kubo tc calculating ------ ##
	vol = read_vol(path+lammpsdata)
	Jx,Jy,Jz = read_heatflux(path+logfile, 'v_Jx','v_Jy','v_Jz',vol,nlog=1)

	plot_hfacf_intersect(hfacf_file,Jx,Jy,Jz,sample_interval,correlate_step,thermointer,timestep)

	green_kubo_tc(hfacf_file,runingtc_file,sample_interval,timestep,TK,vol,ave_factor=0.8)
	
	# ----------- Set the end time and print the execution time -----------
	t_end = time.perf_counter()
	wall_time = convert_time(t_end - t_start)
	print(40*'-')
	print("wall-time: " + str(wall_time))
	print(40*'-')	