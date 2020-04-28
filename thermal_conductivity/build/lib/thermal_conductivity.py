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
	with open(NPT_data,'r')as data,open('log.txt','a')as log:
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
					# print('Size of x direction of system =',system_size_x)
				elif line[2]=='ylo':
					yhi = float(line[1])
					ylo = float(line[0])					
					system_size_y = (yhi-ylo)/10#nm
					# print('Size of y direction of system =',system_size_y)
				# elif line[2]=='zlo':
				# 	system_size_z = (float(line[1])-float(line[0]))/10#nm
				# 	print('Size of z direction of system =',system_size_z)
		size_layer = system_size_x/number_layers
		print('x = ',system_size_x,';','y = ',system_size_y,file=log)
		print('size of everylayer = ',size_layer,file=log)
		print('xhi = ',xhi,';','xlo = ',xlo,file=log)
		print('yhi = ',yhi,';','ylo = ',ylo,file=log)
		# print('\n',file=log)
		print('size() done!')

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
	return  print('temp_grad() done!')

#------------------Calculate the temperature gradient in different way--------------------#
def temp_grad_dTdL(filename,number_layers,number_fixed,number_bath):
	L = system_size_x-system_size_x/number_layers*(number_fixed*2+number_bath*2) #nm
	# wentidu=0.0
	global wentidu
	with open(filename)as wendu,open('log.txt','a')as log:
		for line in wendu:
			L_line=line.strip().split()
			if float(L_line[0])==number_fixed+number_bath+1:
				print('第',L_line[0],'层，高温：',L_line[2],file=log)
				hight_temp=float(L_line[2])#K
			elif float(L_line[0])==number_layers-number_fixed-number_bath and float(L_line[2])>290:
				print('第',L_line[0],'层，低温：',L_line[2],file=log)
				low_temp=float(L_line[2])#K

		wentidu=(hight_temp-low_temp)/L
		print('Temperature gradient (dT/dL)：',wentidu,file=log)
	return  print('temp_grad_dTdL() done!')

#------------------Read input and output energies for calculating heat flux--------------------#
def heat_flux(filename1,filename2,timestep,J2ev):
	with open(filename1) as ener_1,\
		open(filename2,"w") as Q:

		for index, line in enumerate(ener_1):
			for line in ener_1:
				ener_11=line.strip().split()
				# print(ener_11[1],ener_11[2])
				Time=float(ener_11[0])*timestep
				#将能量平均并将ev转换为为J
				Energy=0.5*(float(ener_11[2])-float(ener_11[1]))*J2ev
				# print(Energy)
				#将时间与热流写入文件
				Q.write(str(ener_11[0]))
				Q.write(" ")
				Q.write(str(Time))
				Q.write(" ")
				Q.write(str(Energy))
				Q.write("\n")
		ener_1.close()
		Q.close()
	return print('heat_flux() done!')

#------------------Plot and fitting temperature profile--------------------#
def plot_temp(filename2,layers_fixed,number_fixed,number_bath,i,Plot=True):
	log = open('log.txt','a')

	temp_12=open(filename2,"r")
	# print(temp_12.read())
	# xmin=5#nm
	# xmax=50#nm
	Temp_xmin=5*(number_fixed+number_bath)*size_layer#nm
	Temp_xmax=system_size_x-Temp_xmin#nm	
	x1=list()
	y1=list()#全部
	x2=list()
	y2=list()#拟合
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
	fit = np.polyfit(x2,y2,1)#用1次多项式拟合
	fit_fn = np.poly1d(fit)
	print("Formula of tmperature grafient Fitting:",fit_fn,file=log)#拟合多项式
	print("Slope:" ,fit[0],"(K/nm)","\n"+"Intercept:",fit[1],"(K)",file=log)
	# print('\n',file=log)
	global Temperature_gradient
	Temperature_gradient=fit[0]
	if Plot==True:
		#-plot
		plt.scatter(x1,y1)
		plt.plot(x2,fit_fn(x2),"r-",linewidth=4.0)
		plt.title("Temperature profile")
		plt.xlabel("Distance (nm)")
		plt.ylabel("Temperature (K)")
		plt.savefig(str(i)+"Temperature profile.png")
		plt.show()
		plt.close()
		log.close()

	return print('plot_temp() done!')
#------------------Plot and fitting heat flux--------------------#
def plot_heatflux(filename2,Energy_xmin,Energy_xmax,i,Plot=True):

	with open(filename2) as Q:
		x1=list()
		y1=list()#全部
		x2=list()
		y2=list()#拟合
		for lines in Q:
			Q_t=list(map(eval,lines.split()))
			x1.append(Q_t[1])
			y1.append(Q_t[2])
	#拟合
			if float(Q_t[1])>=Energy_xmin and float(Q_t[1])<=Energy_xmax:
				x2.append(Q_t[1])
				y2.append(Q_t[2])
	log = open('log.txt','a')
	fit = np.polyfit(x2,y2,1)#用1次多项式拟合
	fit_fn1 = np.poly1d(fit)
	print("Formula of Heat flux Fitting:",fit_fn1,file=log)
	print("Heat flux:",fit[0],"(J/ns)",file=log)
	print("Intercept:",fit[1],"(J)\n",file=log)
	# print('\n',file=log)
	global Heat_flux
	Heat_flux=fit[0]
	if Plot == True:
		#-plot
		plt.plot(x1,y1,"o",linewidth=9.0)
		plt.plot(x2,fit_fn1(x2),"r-",linewidth=3.0)
		plt.title("Heat flux (J/ns)")
		plt.xlabel("Time (ns)")
		plt.ylabel("Energy (J)")
		plt.savefig(str(i)+"Heat flux.png")
		plt.show()
		plt.close()
		log.close()
	return Heat_flux, print('plot_heatflux() done!')

#------------------Calculate thermal conductivity--------------------#
def Thermal_conductivity(filename3,thickness,i,dTdL=False):
#Thermal conductivities are saved in filename3
	with open(filename3,"a+") as tc_k,open('log.txt','a')as log:
		A=system_size_y*thickness*(1e-18)#(m2)
		# def TC(Heat_flux,Temperature_gradient,A=2.85e-18):
		if dTdL==True:
			k=Heat_flux/(A*wentidu)
		else:
			k=-Heat_flux/(A*Temperature_gradient)

		print('Thermal conductivity:'+str(round(k,4)),'W/m-K\n',file=log)#round(要输出的值,保留几位小数)
		print('\n',file=log)
		tc_k.write(str(round(k,4)))
		if i == 3:
			tc_k.write('\n')
		else:
			tc_k.write(' ')
	return print('Thermal_conductivity() done!')

#write data for calculating the Transmission
#default MoS2, if number_atom_type==3, it is MoS2-MoSe2 heterostructure
def write_dataforTransmission(filename1,filename2,number_atom_type=2):
	global xhi
	with open(filename1,'r')as former, open(filename2,'w')as latter:
		# data = former.read()
		# print(data)
		for index,line in enumerate(former,1):
			if index<=25:
				latter.write(line)	
			line = line.strip().split()
			# print(line)
			size_line = len(line)
			# print(size_line)
			if size_line==12:
				if float(line[4])<=xhi/2:
					#If the x-coordinate less or equal than half of xhi, 
					#write pristine coordinate
					latter.write('      ')
					latter.write(line[0])
					latter.write('      ')
					latter.write(line[1])
					latter.write('      ')
					latter.write(line[2])
					latter.write('   ')
					latter.write(line[3])
					latter.write('   ')
					latter.write(line[4])
					latter.write('   ')
					latter.write(line[5])
					latter.write('   ')
					latter.write(line[6])
					latter.write('   ')
					latter.write(line[7])
					latter.write('   ')
					latter.write(line[8])
					latter.write('   ')
					latter.write(line[9])
					latter.write('   ')
					latter.write(line[10])
					latter.write(line[11])
					latter.write('\n')
				elif float(line[4])>xhi/2:
					#If it larger than half of xhi, 
					#judge the line[2](type of atom) whether == 1, is S atom
					if number_atom_type==3:
						if line[2]=='1':
							print('S')
							latter.write('      ')
							latter.write(line[0])
							latter.write('      ')
							latter.write(line[1])
							latter.write('      ')
							latter.write(str(4))
							latter.write('   ')
							latter.write(line[3])
							latter.write('   ')
							latter.write(line[4])
							latter.write('   ')
							latter.write(line[5])
							latter.write('   ')
							latter.write(line[6])
							latter.write('   ')
							latter.write(line[7])
							latter.write('   ')
							latter.write(line[8])
							latter.write('   ')
							latter.write(line[9])
							latter.write('   ')
							latter.write(line[10])
							latter.write(line[11])
							latter.write('\n')
						#judge the line[2](type of atom) whether == 2, is Mo atom						
						elif line[2]=='2':
							print('Mo')
							latter.write('      ')
							latter.write(line[0])
							latter.write('      ')
							latter.write(line[1])
							latter.write('      ')
							latter.write(str(5))
							latter.write('   ')
							latter.write(line[3])
							latter.write('   ')
							latter.write(line[4])
							latter.write('   ')
							latter.write(line[5])
							latter.write('   ')
							latter.write(line[6])
							latter.write('   ')
							latter.write(line[7])
							latter.write('   ')
							latter.write(line[8])
							latter.write('   ')
							latter.write(line[9])
							latter.write('   ')
							latter.write(line[10])
							latter.write(line[11])
							latter.write('\n')						
						elif line[2]=='3':
							print('Se')
							latter.write('      ')
							latter.write(line[0])
							latter.write('      ')
							latter.write(line[1])
							latter.write('      ')
							latter.write(str(6))
							latter.write('   ')
							latter.write(line[3])
							latter.write('   ')
							latter.write(line[4])
							latter.write('   ')
							latter.write(line[5])
							latter.write('   ')
							latter.write(line[6])
							latter.write('   ')
							latter.write(line[7])
							latter.write('   ')
							latter.write(line[8])
							latter.write('   ')
							latter.write(line[9])
							latter.write('   ')
							latter.write(line[10])
							latter.write(line[11])
							latter.write('\n')

					elif number_atom_type==2:

						if line[2]=='1':
							print('S')
							latter.write('      ')
							latter.write(line[0])
							latter.write('      ')
							latter.write(line[1])
							latter.write('      ')
							latter.write(str(3))
							latter.write('   ')
							latter.write(line[3])
							latter.write('   ')
							latter.write(line[4])
							latter.write('   ')
							latter.write(line[5])
							latter.write('   ')
							latter.write(line[6])
							latter.write('   ')
							latter.write(line[7])
							latter.write('   ')
							latter.write(line[8])
							latter.write('   ')
							latter.write(line[9])
							latter.write('   ')
							latter.write(line[10])
							latter.write(line[11])
							latter.write('\n')
						#judge the line[2](type of atom) whether == 2, is Mo atom						
						elif line[2]=='2':
							print('Mo')
							latter.write('      ')
							latter.write(line[0])
							latter.write('      ')
							latter.write(line[1])
							latter.write('      ')
							latter.write(str(4))
							latter.write('   ')
							latter.write(line[3])
							latter.write('   ')
							latter.write(line[4])
							latter.write('   ')
							latter.write(line[5])
							latter.write('   ')
							latter.write(line[6])
							latter.write('   ')
							latter.write(line[7])
							latter.write('   ')
							latter.write(line[8])
							latter.write('   ')
							latter.write(line[9])
							latter.write('   ')
							latter.write(line[10])
							latter.write(line[11])
							latter.write('\n')	
	return print('write_dataforTransmission() done!')

#Calculating DeltaT for Transmission
def DeltaT_calculate():
	DeltaT = abs(Temperature_gradient*system_size_x)
	print('DeltaT=',round(DeltaT,4),'K')
	return DeltaT,print('DeltaT_calculate() done!')