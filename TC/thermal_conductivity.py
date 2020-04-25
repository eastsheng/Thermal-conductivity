import matplotlib.pyplot as plt
import numpy as np

#-------The size of system are read from NPT_data in MD simulation-------#
def size(NPT_data,number_layers):
	global system_size_x
	global system_size_y
	global size_layer
	global xlo,xhi
	global ylo,yhi
	global zhi,zhi
	with open(NPT_data,'r')as data:
		for line in data:
			# print(line)
			line = line.strip().split()
			length_line = len(line)
			# print(length_line)
			if length_line==4 and line[2] in ['xlo','ylo','zlo']:
				# print(line)
				if line[2]=='xlo':
					xhi = float(line[1])
					xlo = float(line[0])
					system_size_x = (xhi-xlo)/10#nm
					print('Size of x direction of system =',system_size_x)
				elif line[2]=='ylo':
					yhi = float(line[1])
					ylo = float(line[0])					
					system_size_y = (yhi-ylo)/10#nm
					print('Size of y direction of system =',system_size_y)
				# elif line[2]=='zlo':
				# 	system_size_z = (float(line[1])-float(line[0]))/10#nm
				# 	print('Size of z direction of system =',system_size_z)
	size_layer = system_size_x/number_layers
	print('x = ',system_size_x,';','y = ',system_size_y,'\n')
	print('size of everylayer = ',size_layer)
	print('xhi = ',xhi,';','xlo = ',xlo)
	print('yhi = ',yhi,';','ylo = ',ylo)

	return	system_size_x, system_size_y, size_layer, xhi, xlo, yhi, ylo#,system_size_z



#------------------Read temperature profile for calculating temperature gradient--------------------#
def temp_grad(filename1,filename2,number_layers):

	with open(filename1) as temp_11,\
		open(filename2,"w") as temp_12:
		for index, line in enumerate(temp_11,1):
			temp_gradient = line.strip().split()
			len_temp = len(temp_gradient)
			# print(len_temp)
			if len_temp == 4 and temp_gradient[0] is not '#' :
				coord = float(temp_gradient[0]) #x(nm)
				coord1=round(coord*(system_size_x/number_layers),4)
				temperature = round(float(temp_gradient[3]),4)#T(K)
				# print(coord,temperature)
				temp_12.write(str(coord))
				temp_12.write(" ")
				temp_12.write(str(coord1))
				temp_12.write(" ")
				temp_12.write(str(temperature))
				temp_12.write("\n")
	return 

#------------------Calculate the temperature gradient in different way--------------------#
def temp_grad_dLdT(filename,number_layers,number_fixed,number_bath):
	L = system_size_x-system_size_x/number_layers*(number_fixed*2+number_bath*2) #nm
	# wentidu=0.0
	global wentidu
	with open(filename)as wendu:
		for line in wendu:
			L_line=line.strip().split()
			if float(L_line[0])==number_fixed+number_bath+1:
				print('Number',L_line[0],'layer，high temperature：',L_line[2])
				hight_temp=float(L_line[2])#K
			elif float(L_line[0])==number_layers-number_fixed-number_bath and float(L_line[2])>290:
				print('Number',L_line[0],'layer，low temperature：',L_line[2])
				low_temp=float(L_line[2])#K

		wentidu=(hight_temp-low_temp)/L
		print('Temperature gradient (dL/dT)：',wentidu)
	return

#------------------Read input and output energies for calculating heat flux--------------------#
def heat_flux(filename1,filename2,timestep,J2ev):
	with open(filename1) as ener_1,\
		open(filename2,"w") as Q:

		for index, line in enumerate(ener_1):
			for line in ener_1:
				ener_11=line.strip().split()
				# print(ener_11[1],ener_11[2])
				Time=float(ener_11[0])*timestep
				Energy=0.5*(float(ener_11[2])-float(ener_11[1]))*J2ev
				# print(Energy)
				Q.write(str(ener_11[0]))
				Q.write(" ")
				Q.write(str(Time))
				Q.write(" ")
				Q.write(str(Energy))
				Q.write("\n")
		ener_1.close()
		Q.close()
	return 

#------------------Plot and fitting temperature profile--------------------#
def plot_temp(filename2,layers_fixed,number_fixed,number_bath,i):

	temp_12=open(filename2,"r")
	Temp_xmin=4*(number_fixed+number_bath)*size_layer#nm
	Temp_xmax=system_size_x-Temp_xmin#nm	
	x1=list()
	y1=list()
	x2=list()
	y2=list()
	for lines in temp_12:
		temp_gradient1 = lines.strip().split()
		temp_gradient  = list(map(eval,temp_gradient1))
		# print(temp_gradient)
		if temp_gradient[0] not in layers_fixed:

			x1.append(temp_gradient[1])
			y1.append(temp_gradient[2])

			if float(temp_gradient[1])>=Temp_xmin and float(temp_gradient[1])<=Temp_xmax:
				x2.append(temp_gradient[1])
				y2.append(temp_gradient[2])
				# print(type(temp_gradient[1]))
	fit = np.polyfit(x2,y2,1)
	fit_fn = np.poly1d(fit)
	print("Fitting:",fit_fn)
	print("slope-temperature gradient:" ,fit[0],"(K/nm)","\n"+"intercept:",fit[1],"(K)","\n")
	global Temperature_gradient
	Temperature_gradient=fit[0]
	#-------plot
	plt.scatter(x1,y1)
	plt.plot(x2,fit_fn(x2),"r-",linewidth=4.0)
	plt.title("Temperature profile")
	plt.xlabel("Distance (nm)")
	plt.ylabel("Temperature (K)")
	plt.savefig(str(i)+"Temperature profile.png")
	plt.show()
	plt.close()
	return
#------------------Plot and fitting heat flux--------------------#
def plot_heatflux(filename2,Energy_xmin,Energy_xmax,i):

	with open(filename2) as Q:
		x1=list()
		y1=list()
		x2=list()
		y2=list()
		for lines in Q:
			Q_t=list(map(eval,lines.split()))
			x1.append(Q_t[1])
			y1.append(Q_t[2])
	#Fitting
			if float(Q_t[1])>=Energy_xmin and float(Q_t[1])<=Energy_xmax:
				x2.append(Q_t[1])
				y2.append(Q_t[2])
	fit = np.polyfit(x2,y2,1)
	fit_fn1 = np.poly1d(fit)
	print("Fitting：",fit_fn1)
	print("Heat flux：",fit[0],"(J/ns)")
	print("Intercept：",fit[1],"(J)\n")
	global Heat_flux
	Heat_flux=fit[0]
	#-------plot
	plt.plot(x1,y1,"o",linewidth=9.0)
	plt.plot(x2,fit_fn1(x2),"r-",linewidth=3.0)
	plt.title("Heat flux (J/ns)")
	plt.xlabel("Time (ns)")
	plt.ylabel("Energy (J)")
	plt.savefig(str(i)+"Heat flux.png")
	plt.show()
	plt.close()
	return Heat_flux

#------------------Calculate thermal conductivity--------------------#
def Thermal_conductivity(filename3,thickness,i,dLdT=False):
#Thermal conductivities are saved in filename3
	with open(filename3,"a+") as tc_k:
		A=system_size_y*thickness*(1e-18)#(m2)
		# def TC(Heat_flux,Temperature_gradient,A=2.85e-18):
		if dLdT==True:
			k=Heat_flux/(A*wentidu)
		else:
			k=-Heat_flux/(A*Temperature_gradient)

		print('thermal conductivity:'+str(round(k,4)),'W/m-K\n')
		tc_k.write(str(round(k,4)))
		if i == 3:
			tc_k.write('\n')
		else:
			tc_k.write(' ')
	return


