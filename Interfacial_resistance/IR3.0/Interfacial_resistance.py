#For calcualting Interfacial resistance
from scipy import integrate
import numpy as np 
import matplotlib.pyplot as plt
import fastdataing as fd
import readlammpsdata as rld
from print_log import print_log

def Interfacial_resistance(TEfile,fit_ranges=[400,2],area=False,timestep=1):
	'''
	Interfacial resistance is calculated by integral change of temperature difference and recording total energy.
	Parameters:
	- TEfile: Temperature and energy file, # TimeStep c_TH c_TQ c_allke c_allpe
	- fit_ranges: default = fit_ranges=[400,2]
	- area: cross sectional area
	- timestep: fs
	'''
	kcalmol2J = 4186.0/6.02214e23
	ev2J = 1.60217662e-19 # ev2J
	ps2s,fs2s = 1e-12, 1e-15
	energy2J = ev2J 
	time2s = ps2s 

	data = np.loadtxt(TEfile)
	inter_step = int(data[1,0] - data[0,0])
	time = (data[:,0]-data[0,0])*timestep*time2s # s
	inter_time = inter_step*timestep*time2s
	step = np.array([i for i in range(len(data))])
	Temperature_up, Temperature_down = data[:,1], data[:,2]
	ke, pe = data[:,3], data[:,4]
	total_energy = (ke+pe)*energy2J

	#**********define temperature difference function**********#
	def Temp_Diff(x):
		temp_diff=(Temperature_up[x]-Temperature_down[x])
		return temp_diff
	#**********Integrate**********#	
	y= Temp_Diff(step)
	DT_integrate= integrate.cumulative_trapezoid(y)

	x1 = DT_integrate
	y1 = total_energy[1:]
	#******control the fitting interval******#
	Fit_minx2 = int(len(step)/fit_ranges[0])
	Fit_maxx2 = int(len(step)/fit_ranges[1])
	x2 = DT_integrate[Fit_minx2:Fit_maxx2]
	y2 = total_energy[Fit_minx2:Fit_maxx2]

	fit = np.polyfit(x2,y2,1)
	fit_fn = np.poly1d(fit)
	print("Formula of Heat flux Fitting:y2 = ",fit_fn)
	#Fitting slope
	Area_R = fit[0]/inter_time # change unit to K.s
	R = -area/Area_R
	print(f'R = {R} Km^2/W\n')
	return  time, Temperature_up, Temperature_down, total_energy, x1,y1,x2,fit_fn(x2)


if __name__ == '__main__':
	
	fit_ranges=[500,5]#不用怎么变
	path = "./"
	logfile = f'{path}calc_log.log'
	print_log(logfile)
	k = 1
	for i in range(1,k+1):
		lmp = f'{path}GRA_C3N_npt1.data'
		TEfile = f'{path}A3relaxation'+str(i)+'.dat'
		x = rld.read_len(lmp,"x")
		y = rld.read_len(lmp,"y")
		area = x*y*1e-20 # m^2
		print(f">>> Area = {area} m^2")
		time, Temperature_up, \
		Temperature_down, total_energy, \
		x1,y1,x2,y2 \
			= Interfacial_resistance(
			TEfile=TEfile,fit_ranges=fit_ranges,area=area,timestep=5e-4)


	fig = fd.add_fig(figsize=(17.5,6))
	plt.subplots_adjust(wspace=0.45)
	ax = fd.add_ax(fig,subplot=(1,2,1))
	ax.plot(time,Temperature_up,c="r",alpha=0.6,label="Temperature 1")
	ax.plot(time,Temperature_down,c="g",alpha=0.6,label="Temperature 2")
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Temperature (K)")
	fd.set_fig(ax,loc="upper left")	
	ay = ax.twinx()
	ay.plot(time,total_energy,c="b",alpha=0.6,label="Total energy")
	ay.set_ylabel("Energy (J)")
	fd.set_fig(ay)
	ax = fd.add_ax(fig,subplot=(1,2,2))
	ax.plot(x1,y1,c="gray",alpha=0.8,linewidth=6.0)
	ax.plot(x2,y2,"r-",linewidth=3.0)
	ax.set_xlabel("DT (K.step)")
	ax.set_ylabel("Energy (J)")
	plt.savefig("IR.png",dpi=300)
	plt.show()
