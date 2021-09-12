#For calcualting Interfacial resistance
from scipy import integrate
import numpy as np 
import matplotlib.pyplot as plt
#-------The size of system are read from NPT_data in MD simulation-------#
def Area(NPT_data,logfile):
	'''Area of out of plane is obtained from NPT data.'''
	with open(NPT_data,'r')as data,open(logfile,'a')as log:
		log.write('********Run'+str(i)+'********\n')
		for line in data:
			line = line.strip().split()
			length_line = len(line)
			# print(length_line)
			if length_line==4 and line[2] in ['xlo','ylo','zlo']:
				# print(line)
				if line[2]=='xlo':
					xhi = float(line[1])
					xlo = float(line[0])
					system_size_x = (xhi-xlo)/10#nm
					print('Size of x direction of system =',system_size_x,file=log)
				elif line[2]=='ylo':
					yhi = float(line[1])
					ylo = float(line[0])					
					system_size_y = (yhi-ylo)/10#nm
					print('Size of y direction of system =',system_size_y,file=log)
				# elif line[2] == 'zlo':
				# 	zhi = float(line[1])
				# 	zlo = float(line[0])
				# 	system_size_z = (yhi-ylo)/10#nm
				# 	print('Size of z direction of system =',system_size_z,file=log)
		global area
		area = system_size_x*system_size_y*(1e-18)
		print('Area',str(i),' = ',area,'m^2',file=log)
	return	 print('\n**********Area done!**********\nArea',str(i),'= ',area,'m^2\n')

def Interfacial_resistance(Temp_and_energy_file,logfile,IR_file,\
	fit_range1=120,fit_range2=4,timestep=5e-4,inter_step=100,Figure=True):
	'''Interfacial resistance is calculated by integral change of temperature difference 
	and recording total energy.
	Format of each line:
	step Temperature_up Temperature_down kinetic_energy potential_energy'''
	#**********variable**********#
	global i 
	inter_time = timestep*inter_step*(1e-12)#s
	ev2J = 1.60217662e-19#ev2J
	#**********read data**********#
	data = np.loadtxt(Temp_and_energy_file)
	# print(data.shape)
	log = open(logfile,'a')
	#**********output of result**********#
	IR = open(IR_file,'a')#Interfacial Resistance
	#**********step 2 time**********#
	step = (data[:,0]-data[0,0])/inter_step
	step = list(map(int,step))
	step = np.array(step)
	time = inter_time*step#s
	#**********Plot temperature profile**********#
	Temperature_up = data[:,1]#high temperature
	Temperature_down = data[:,2]#low temperature
	plt.figure(num=1,figsize=(8,6))
	plt.plot(time,Temperature_up,time,Temperature_down)
	plt.title("Temperature")
	plt.xlabel("Time (s)")
	plt.ylabel("Temperature (K)")
	plt.savefig(str(i)+"Temperature profile.png")
	if Figure == True:
		plt.show()
	plt.close()
	#**********Plot Total energy profile**********#
	kinetic_energy = data[:,3]
	potential_energy = data[:,4]
	total_energy = (kinetic_energy+potential_energy)*ev2J
	plt.figure(num=2,figsize=(8,6))	
	plt.plot(time,total_energy)
	plt.title("Total energy")
	plt.xlabel("Time (s)")
	plt.ylabel("Energy (J)")
	plt.savefig(str(i)+"Total energy profile.png")	
	if Figure==True:
		plt.show()
	plt.close()
	# print(step)

	#**********define temperature difference function**********#
	def Temp_Diff(x):
		temp_diff=(Temperature_up[x]-Temperature_down[x])
		return temp_diff
	#**********Integrate**********#	
	y= Temp_Diff(step)
	DT_integrate= integrate.cumtrapz(y)

	#**********Fitting and plot**********#	
	x1 = DT_integrate
	y1 = total_energy[1:]
	#******control the fitting interval******#
	Fit_minx2 = int(len(step)/fit_range1)
	# print(Fit_minx2)
	Fit_maxx2 = int(len(step)/fit_range2)
	# print(Fit_maxx2)
	x2 = DT_integrate[Fit_minx2:Fit_maxx2]
	y2 = total_energy[Fit_minx2:Fit_maxx2]
	# x2 = DT_integrate[20:1600]
	# y2 = total_energy[20:1600]

	fit = np.polyfit(x2,y2,1)
	fit_fn1 = np.poly1d(fit)
	print("Formula of Heat flux Fitting:y2 = ",fit_fn1,file=log)
	#Fitting slope
	Area_R = fit[0]/inter_time#change unit to K.s
	R = -area/Area_R
	print('R'+str(i)+'=',R,'Km^2/W\n',file=log)

	plt.figure(num=3,figsize=(8,6))
	plt.plot(x1,y1,linewidth=6.0)
	plt.plot(x2,fit_fn1(x2),"r-",linewidth=3.0)
	plt.title("Total energy")
	plt.xlabel("DT (K.step)")
	plt.ylabel("Energy (J)")
	plt.savefig(str(i)+"Total energy profile-DT.png")	
	if Figure==True:
		plt.show()
	plt.close()
	#output Interfacial resistance
	IR.write(str(R))
	IR.write('  ')
	IR.close()
	log.close()
	return print('**********Interfacial_resistance done!**********'),\
	print('R'+str(i)+'=',R,'Km^2/W\n\n')


#**********Main Program**********#
fit_range1=500
fit_range2=4

k = 3
for i in range(1,k+1):
#Area(nptfile,logfile)
	Area('Bulk_MoS2_npt'+str(i)+'.data','log.txt')
#Interfacial resistance(temperature and energy file, logfile, result file,
#timestep,interval of step,whether show figure)
	Interfacial_resistance('A3relaxation'+str(i)+'.dat',\
		'log.txt',\
		'Interfacial_resistance.result',\
		fit_range1,fit_range2,timestep=5e-4,inter_step=100,Figure=True)


