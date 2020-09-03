#Calculating the thermal condcutivity of 2D materials from NEMD dump------updated version2.0
import matplotlib.pyplot as plt
import numpy as np
'''
if heatflux_direction=1,namely, x direction is the heat flux direction
elif heatflux_direction=2,namely, y direction is the heat flux direction
'''
class ThermalConductivity(object):

	def read_size(self,relax_data,i,thickness,heatflux_direction=1):
		self.relax_data =relax_data
		self.case = i
		self.system_size_x = 0 #初始化
		self.system_size_y = 0
		self.system_size_z = 0
		self.thickness = thickness
		self.area = 0 #初始化
		self.heatflux_direction = heatflux_direction

		with open(self.relax_data,'r')as data:
			for line in data:
				line = line.strip().split()
				length_line = len(line)
				if length_line==4 and line[2] in ['xlo','ylo','zlo']:
					if line[2]=='xlo':
						xhi = float(line[1])
						xlo = float(line[0])
						self.system_size_x = (xhi-xlo)/10
					elif line[2]=='ylo':
						yhi = float(line[1])
						ylo = float(line[0])						
						self.system_size_y = (yhi-ylo)/10
					else:
						line[2] == 'zlo'
						zhi = float(line[1])
						zlo = float(line[0])
						self.system_size_z = (zhi-zlo)/10#nm 真空层，一般用不到
			if heatflux_direction==1:# 如果热流方向为x，y即是宽度
				self.area = self.system_size_y*self.thickness*(1e-18)
			elif heatflux_direction==2:
				self.area = self.system_size_x*self.thickness*(1e-18)	
			print('Area_'+str(self.case),'=',self.area,'m^2\n')
		return
		
	#------------------Read temperature profile for calculating temperature gradient--------------------#
	def temp_grad(self,tempfile,number_layers,number_fixed,number_bath,fit_factor,Plot=True):
		self.tempfile = tempfile
		self.number_layers = number_layers
		self.number_fixed = number_fixed
		self.number_bath = number_bath
		self.fit_factor = fit_factor
		self.plot = Plot
		self.Temperature_gradient_fit = 0 #初始化,拟合所得温度梯度
		self.Temperature_gradient_difference1 = 0 #初始化,温差所得温度梯度（不包括热浴层的温度）
		self.Temperature_gradient_difference2 = 0 #初始化,直接初始温差所得温度梯度

		temp_data=np.loadtxt(self.tempfile,skiprows=4)
		if self.heatflux_direction==1:
			thickness_eachlayer = self.system_size_x/self.number_layers
		elif self.heatflux_direction==2:
			thickness_eachlayer = self.system_size_y/self.number_layers

		coord_x=temp_data[:,0]*(thickness_eachlayer)#根据自己的数据修改
		temperature=temp_data[:,3]
		# plot x1,y1
		x1=coord_x[self.number_fixed:self.number_layers-self.number_fixed]#不包括固定层
		y1=temperature[self.number_fixed:self.number_layers-self.number_fixed]
		# 拟合范围
		fit_range_T = (self.number_fixed+self.number_bath)*self.fit_factor
		# 拟合的x2,y2
		x2=coord_x[fit_range_T:self.number_layers-fit_range_T]
		y2=temperature[fit_range_T:self.number_layers-fit_range_T]
		fit=np.polyfit(x2,y2,1)
		fit_fn = np.poly1d(fit)
		# 拟合的温度梯度
		self.Temperature_gradient_fit=fit[0]	

		print("Formula of tmperature grafient Fitting: y= ",fit_fn)#拟合多项式
		print("Slope:" ,self.Temperature_gradient_fit,"(K/nm)","\n"+"Intercept:",fit[1],"(K)")
		# plot
		plt.rc('font', family='Times New Roman', size=16)
		plt.scatter(x1,y1)
		plt.plot(x2,fit_fn(x2),'r-',linewidth=4.0)
		plt.title("Temperature profile")
		plt.xlabel("x coord (nm)")
		plt.ylabel("Temperature (K)")
		plt.savefig(str(self.case)+"Temperature profile.png")
		if Plot==True:
			plt.show()
		plt.close()
		# 第二种温度梯度,除去固定层和热浴层
		L1=self.system_size_x-thickness_eachlayer*(self.number_fixed*2+self.number_bath*2)
		high_temp1=temperature[self.number_fixed+self.number_bath]
		low_temp1=temperature[self.number_layers-self.number_fixed-self.number_bath-1]
		self.Temperature_gradient_difference1=(high_temp1-low_temp1)/L1
		# 第二种温度梯度，直接使用初始温差,包括热浴层
		L2=self.system_size_x-thickness_eachlayer*(self.number_fixed*2)
		high_temp2=temperature[self.number_fixed]
		low_temp2=temperature[self.number_layers-self.number_fixed-1]
		self.Temperature_gradient_difference2=(high_temp2-low_temp2)/L2 

		print('\nTemperature_gradient_difference1 = ',self.Temperature_gradient_difference1)
		print('High low temperature1(K)',high_temp1,low_temp1)
		print('size L1',L1)

		print('\nTemperature_gradient_difference2 = ',self.Temperature_gradient_difference2)
		print('High low temperature2(K)',high_temp2,low_temp2)
		print('size L2',L2)

		return 


	#------------------Read input and output energies for calculating heat flux--------------------#
	def heat_flux(self,energyfile,timestep=5e-7):
		self.energyfile = energyfile
		self.timestep = timestep#ns
		self.J2ev = 1.602763e-19
		self.Heat_flux = 0

		energy_data = np.loadtxt(self.energyfile,skiprows=2)

		time = energy_data[:,0]*self.timestep
		energy = 0.5*(energy_data[:,2]-energy_data[:,1])*self.J2ev

		x1 = time
		y1 = energy

		fit_range_E = len(x1)

		x2 = time[2:fit_range_E-2]
		y2 = energy[2:fit_range_E-2]

		fit = np.polyfit(x2,y2,1)
		fit_fn = np.poly1d(fit)

		print("\nFormula of Heat flux Fitting: y = ",fit_fn)
		print("Heat flux:",fit[0],"(J/ns)")
		print("Intercept:",fit[1],"(J)\n")

		self.Heat_flux = fit[0]
		plt.rc('font', family='Times New Roman', size=16)
		plt.plot(x1,y1,"o",linewidth=9.0)
		plt.plot(x2,fit_fn(x2),"r-",linewidth=3.0)
		plt.title("Heat flux (J/ns)")
		plt.xlabel("Time (ns)")
		plt.ylabel("Energy (J)")
		plt.savefig(str(self.case)+"Heat flux.png")
		if self.plot==True:
			plt.show()
		plt.close()
		return 

	'''   
	TempGrad_fator=1,use fitting temperature gradient.
	TempGrad_fator=2,without including highest and lowest temperatures,namely hot and cold bath.
	TempGrad_fator=3,use directly temperature difference.
	'''
	def thermal_conductivity(self,result,TempGrad_fator=1):
		self.result = result
		self.TempGrad_fator = TempGrad_fator
		self.k = 0

		with open(self.result,"a+") as tc_k:
			if TempGrad_fator==1:
				self.k=abs(self.Heat_flux/(self.area*self.Temperature_gradient_fit))
			elif TempGrad_fator==2:
				self.k=abs(self.Heat_flux/(self.area*self.Temperature_gradient_difference1))
			elif TempGrad_fator==3:
				self.k=abs(self.Heat_flux/(self.area*self.Temperature_gradient_difference2))
			else:
				print('\n********TempGrad_fator is wrong!********\n')

			print('Thermal conductivity:'+str(round(self.k,4)),'W/m-K\n')#round(要输出的值,保留几位小数)
			tc_k.write(str(round(self.k,4)))
			if self.case == 3:
				tc_k.write('\n')
			else:
				tc_k.write(' ')
		print('\n**********Thermal Conductivity Calculations are Completed**********\n')		
		return 	

	def logfile(self,logname):
		self.logname = logname
		with open(logname,'a') as log:
			log.write(5*'-'+' Run '+str(self.case)+5*'-'+'\n')

			log.write('x = '+str(self.system_size_x)+' nm\n')
			log.write('y = '+str(self.system_size_y)+' nm\n')
			log.write('z = '+str(self.system_size_z)+' nm\n')
			log.write('Area = '+str(self.area)+' m^2\n')

			log.write('Temperature gradient of fit = '+str(self.Temperature_gradient_fit))
			log.write('\nTemperature gradient of difference1 = '+str(self.Temperature_gradient_difference1))
			log.write('\nTemperature gradient of difference2 = '+str(self.Temperature_gradient_difference2))

			log.write('\nHeat Flux = '+str(self.Heat_flux)+'(J/ns)')

			log.write('\nThermal Conductivity = '+str(self.k)+'(W/m-K)\n')
			

# i = 1
# System_temp=300
# number_layers = 120
# number_fixed = 2
# number_bath = 2
# heatflux_direction = 1
# thickness=0.61
# timestep = 5e-7#ns

# TC = ThermalConductivity()

# TC.read_size('./0811/TC3.0/MoS2_NPT.data',i,thickness,heatflux_direction)

# temperaturefile = './0811/TC3.0/'+str(i)+"_temp_equ_"+str(System_temp)+"K.dat"
# TC.temp_grad(temperaturefile,number_layers,number_fixed,number_bath,fit_factor=3,Plot=True)


# heatfluxfile = './0811/TC3.0/'+str(i)+"_Ener_equ_"+str(System_temp)+"K.dat"
# TC.heat_flux(heatfluxfile,timestep)

# tcfile = "./0811/TC3.0/Thermal_conductivity.txt"
# TC.thermal_conductivity(tcfile,TempGrad_fator=1)


# TC.logfile('./0811/TC3.0/log.txt')