#For calcualting Interfacial resistance
from scipy import integrate
import numpy as np 
import matplotlib.pyplot as plt
class InterfacialResistance(object):
	#-------The size of system are read from relax_data in MD simulation-------#
	def read_size(self, relax_data, i):
		self.relax_data = relax_data
		self.case = i
		self.system_size_x = 0  # 初始化
		self.system_size_y = 0
		self.system_size_z = 0
		self.area = 0  # 初始化
		with open(self.relax_data, 'r')as data:
			for line in data:
				line = line.strip().split()
				length_line = len(line)
				if length_line == 4 and line[2] in ['xlo', 'ylo', 'zlo']:
					if line[2] == 'xlo':
						xhi = float(line[1])
						xlo = float(line[0])
						self.system_size_x = (xhi-xlo)/10
					elif line[2] == 'ylo':
						yhi = float(line[1])
						ylo = float(line[0])
						self.system_size_y = (yhi-ylo)/10
					else:
						line[2] == 'zlo'
						zhi = float(line[1])
						zlo = float(line[0])
						self.system_size_z = (zhi-zlo)/10 

			self.area = self.system_size_x*self.system_size_y*(1e-18)
			print('Area_'+str(self.case), '=', self.area, 'm^2\n')
		return

	def Interfacial_resistance(self,Temp_and_energy_file,IR_file,savepath,
                            fit_range1, fit_range2, timestep=5e-4, inter_step=100, Figure=True):
		self.Temp_and_energy_file = Temp_and_energy_file
		self.IR_file = IR_file
		self.fit_range1 = fit_range1
		self.fit_range2 = fit_range2
		self.timestep = timestep  # ps
		self.inter_step = inter_step # step
		self.ev2J = 1.60217662e-19  # ev2J
		self.Figure = Figure # whether plot
		self.R = 0#初始化热阻
		'''Interfacial resistance is calculated by integral change of temperature difference 
		and recording total energy.
		Format of each line of md dump file:
		step Temperature_up Temperature_down kinetic_energy potential_energy'''
		#**********variable**********#
		inter_time = self.timestep*self.inter_step*(1e-12)#s
		#**********read data**********#
		data = np.loadtxt(self.Temp_and_energy_file)
		# print(data.shape)
		#**********output of result**********#
		IR = open(self.IR_file,'a')#Interfacial Resistance
		#**********step 2 time**********#
		step = (data[:,0]-data[0,0])/self.inter_step # 第一步 步数归零
		step = list(map(int,step))
		step = np.array(step)
		time = inter_time*step#总时间（s）
		#**********Plot temperature profile**********#
		Temperature_up = data[:,1]#high temperature
		Temperature_down = data[:,2]#low temperature
		plt.rc('font', family='Times New Roman', size=16)
		plt.figure(num=1,figsize=(8,6))
		plt.plot(time,Temperature_up,time,Temperature_down)
		plt.title("Temperature")
		plt.xlabel("Time (s)")
		plt.ylabel("Temperature (K)")
		plt.savefig(savepath+str(i)+"Temperature profile.png")
		if self.Figure == True:
			plt.show()
		plt.close()
		#**********Plot Total energy profile**********#
		kinetic_energy = data[:,3]
		potential_energy = data[:,4]
		total_energy = (kinetic_energy+potential_energy)*self.ev2J

		plt.figure(num=2,figsize=(8,6))	
		plt.plot(time,total_energy)
		plt.title("Total energy")
		plt.xlabel("Time (s)")
		plt.ylabel("Energy (J)")
		plt.savefig(savepath+str(i)+"Total energy profile.png")
		if self.Figure==True:
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
		Fit_minx2 = int(len(step)/self.fit_range1)
		# print(Fit_minx2)
		Fit_maxx2 = int(len(step)/self.fit_range2)
		# print(Fit_maxx2)
		x2 = DT_integrate[Fit_minx2:Fit_maxx2]
		y2 = total_energy[Fit_minx2:Fit_maxx2]
		# x2 = DT_integrate[20:1600]
		# y2 = total_energy[20:1600]

		fit = np.polyfit(x2,y2,1)
		fit_fn1 = np.poly1d(fit)
		print("Formula of Heat flux Fitting:y2 = ",fit_fn1)
		#Fitting slope
		Area_R = fit[0]/inter_time#change unit to K.s
		self.R = -self.area/Area_R
		print('R'+str(i)+'=', self.R, 'Km^2/W\n')

		plt.figure(num=3,figsize=(8,6))
		plt.plot(x1,y1,linewidth=6.0)
		plt.plot(x2,fit_fn1(x2),"r-",linewidth=3.0)
		plt.title("Total energy")
		plt.xlabel("DT (K.step)")
		plt.ylabel("Energy (J)")
		plt.savefig(savepath+str(i)+"Total energy profile-DT.png")	
		if self.Figure==True:
			plt.show()
		plt.close()
		#output Interfacial resistance
		IR.write(str(self.R))
		IR.write('  ')
		IR.close()
		print('**********Interfacial_resistance done!**********')
		return 

	def logfile(self, logname):
		self.logname = logname
		with open(logname, 'a') as log:
			log.write(5*'-'+' Run '+str(self.case)+5*'-'+'\n')

			log.write('x = '+str(self.system_size_x)+' nm\n')
			log.write('y = '+str(self.system_size_y)+' nm\n')
			log.write('z = '+str(self.system_size_z)+' nm\n')
			log.write('Area = '+str(self.area)+' m^2\n')

			log.write('IR = '+str(self.R)+' Km^2/W\n')
		return 

#**********Main Program**********#
# case
k = 1
# 步长与输出间隔
timestep = 5e-4
inter_step = 100

# 拟合范围，不是线性的
fit_range1 = 500
fit_range2 = 4

InterResist = InterfacialResistance()
# 保存文件的路径
# savepath = './Interfacial_resistance/IR2.0/' #editor : vscode
savepath = '' #editor : sublime

for i in range(1,k+1):
	InterResist.read_size(savepath+'GRA_C3N_npt'+str(i)+'.data', i)

	InterResist.Interfacial_resistance(savepath+'A3relaxation'+str(i)+'.dat',
                                    savepath+'Interfacial_resistance.txt', 
                                    savepath, fit_range1, fit_range2, timestep,
									inter_step, Figure=False)
	InterResist.logfile(savepath+'IRlog.txt')


