#Calculating the thermal condcutivity of 2D materials from NEMD dump------updated version2.0
import matplotlib.pyplot as plt
import numpy as np
'''
if heatflux_direction=1,namely, x direction is the heat flux direction
elif heatflux_direction=2,namely, y direction is the heat flux direction
'''
def Area(NPT_data,logfile,i,thickness=0.61,heatflux_direction=1):
	'''Area of out of plane is obtained from NPT data.'''
	global area
	global system_size_x
	global system_size_y
	with open(NPT_data,'r')as data,open(logfile,'w')as log:
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
		if heatflux_direction==1:
			area = system_size_y*thickness*(1e-18)
		elif heatflux_direction==2:
			area = system_size_x*thickness*(1e-18)
		
		print('Area',str(i),' = ',area,'m^2',file=log)
	return	 print('\n**********Area Done!**********\nArea',str(i),'= ',area,'m^2\n')


#------------------Read temperature profile for calculating temperature gradient--------------------#
def temp_grad(tempfile,number_layers,number_fixed,number_bath,i,heatflux_direction=1,fit_factor=2,Plot=True):

	log = open("log.txt","w")
	temp_data=np.loadtxt(tempfile,skiprows=4)

	if heatflux_direction==1:
		thickness_eachlayer = system_size_x/number_layers
	elif heatflux_direction==2:
		thickness_eachlayer = system_size_y/number_layers

	coord_x=temp_data[:,0]*(thickness_eachlayer)
	temperature=temp_data[:,3]

	x1=coord_x[number_fixed:number_layers-number_fixed]
	y1=temperature[number_fixed:number_layers-number_fixed]

	fit_range_T = (number_fixed+number_bath)*fit_factor

	x2=coord_x[fit_range_T:number_layers-fit_range_T]
	y2=temperature[fit_range_T:number_layers-fit_range_T]
	fit=np.polyfit(x2,y2,1)
	fit_fn = np.poly1d(fit)

	global Temperature_gradient_fit
	Temperature_gradient_fit=fit[0]	

	print("Formula of tmperature grafient Fitting:",fit_fn,file=log)#拟合多项式
	print("Slope:" ,Temperature_gradient_fit,"(K/nm)","\n"+"Intercept:",fit[1],"(K)",file=log)

	plt.scatter(x1,y1)
	plt.plot(x2,fit_fn(x2),'r-',linewidth=4.0)
	plt.title("Temperature profile")
	plt.xlabel("x coord (nm)")
	plt.ylabel("Temperature (K)")
	plt.savefig(str(i)+"Temperature profile.png")
	if Plot==True:
		plt.show()
	plt.close()
	
	global Temperature_gradient_difference1
	global Temperature_gradient_difference2
	
	L1=system_size_x-thickness_eachlayer*(number_fixed*2+number_bath*2)
	high_temp1=temperature[number_fixed+number_bath]
	low_temp1=temperature[number_layers-number_fixed-number_bath-1]
	Temperature_gradient_difference1=(high_temp1-low_temp1)/L1

	L2=system_size_x-thickness_eachlayer*(number_fixed*2)
	high_temp2=temperature[number_fixed]
	low_temp2=temperature[number_layers-number_fixed-1]
	Temperature_gradient_difference2=(high_temp2-low_temp2)/L2   

	print('Temperature_gradient_difference1',high_temp1,low_temp1,file=log)
	print('Temperature_gradient_difference2',high_temp2,low_temp2,file=log)
	print('L1',L1,file=log)
	print('L2',L2,file=log)
	log.close()
	
	return  print('**********Temperature Gradient Done!**********',\
		'\nTemperature_gradient_fit:',Temperature_gradient_fit,\
		'\nTemperature_gradient_difference1:',Temperature_gradient_difference1,\
		'\nTemperature_gradient_difference2:',Temperature_gradient_difference2)


#------------------Read input and output energies for calculating heat flux--------------------#
def heat_flux(energyfile,i,timestep,J2ev=1.602763e-19,Plot=True):

	log = open('log.txt','a')
	energy_data = np.loadtxt(energyfile,skiprows=2)

	time = energy_data[:,0]*timestep
	energy = 0.5*(energy_data[:,2]-energy_data[:,1])*J2ev

	x1 = time
	y1 = energy

	fit_range_E = len(x1)

	x2 = time[2:fit_range_E-2]
	y2 = energy[2:fit_range_E-2]

	fit = np.polyfit(x2,y2,1)
	fit_fn = np.poly1d(fit)


	print("Formula of Heat flux Fitting:",fit_fn,file=log)
	print("Heat flux:",fit[0],"(J/ns)",file=log)
	print("Intercept:",fit[1],"(J)\n",file=log)

	global Heat_flux
	Heat_flux = fit[0]

	plt.plot(x1,y1,"o",linewidth=9.0)
	plt.plot(x2,fit_fn(x2),"r-",linewidth=3.0)
	plt.title("Heat flux (J/ns)")
	plt.xlabel("Time (ns)")
	plt.ylabel("Energy (J)")
	plt.savefig(str(i)+"Heat flux.png")
	if Plot==True:
		plt.show()
	plt.close()
	log.close()

	return print('\n**********Heat Flux Done!**********',\
		"\nHeat flux:",fit[0],"(J/ns)")

'''   
TempGrad_fator=1,use fitting temperature gradient.
TempGrad_fator=2,without including highest and lowest temperatures,namely hot and cold bath.
TempGrad_fator=3,use directly temperature difference.
'''
def Thermal_conductivity(result,i,TempGrad_fator=1):
	with open(result,"a+") as tc_k,open('log.txt','a')as log:
		if TempGrad_fator==1:
			k=-Heat_flux/(area*Temperature_gradient_fit)
		elif TempGrad_fator==2:
			k=Heat_flux/(area*Temperature_gradient_difference1)
		elif TempGrad_fator==3:
			k=Heat_flux/(area*Temperature_gradient_difference2)
		else:
			print('\n********TempGrad_fator is wrong!********\n')

		print('Thermal conductivity:'+str(round(k,4)),'W/m-K\n',file=log)#round(要输出的值,保留几位小数)
		print('\n',file=log)
		tc_k.write(str(round(k,4)))
		if i == 3:
			tc_k.write('\n')
		else:
			tc_k.write(' ')
	return print('\n**********Thermal Conductivity Calculations are Completed**********\n')