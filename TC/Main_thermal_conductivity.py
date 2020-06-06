import thermal_conductivity as TC
##############################################
#case
k = 1
#系统尺寸(nm)
thickness = 0.6469#0.61
#系统温度(K)
System_temp = 300
#层数
number_layers = 120
#固定层数
number_fixed = 2
#热浴层数
number_bath = 2
#固定层
layers_fixed = [1,2,number_layers-1,number_layers,number_layers+1]
#单位换算
timestep=5e-7#ns
J2ev=1.602763e-19#ev转换为J

#热流拟合区间(ns)
Energy_xmin=0.1#ns
Energy_xmax=9.9#ns

#---------------------计算热导率从此处开始计算----------------------#
for i in range(1,k+1):
	#The size of system are read from NPT_data
	#size(npt.data,number of layers)
	TC.size('MoS2_NPT.data',number_layers)

	#Read temperature profile for calculating temperature gradient
	#temp_grad(temperature profile data1,temperature profile data2, number of layers)
	TC.temp_grad(str(i)+"_temp_equ_"+str(System_temp)+"K.dat",\
		str(i)+"_temp_equ_"+str(System_temp)+"K.txt",number_layers)
	#Plot and fitting temperature profile
	#plot_temp(temperature profile data2,layers of fixed,number of fixed,number of bath,i)
	TC.plot_temp(str(i)+"_temp_equ_"+str(System_temp)+"K.txt",\
		layers_fixed,number_fixed,number_bath,i,Templayer_times=5,Plot=True)

	#Read input and output energies for calculating heat flux
	#heat_flux(energy data1,energy data2,timestep,J2ev)
	TC.heat_flux(str(i)+"_Ener_equ_"+str(System_temp)+"K.dat",\
		str(i)+"_Ener_equ_"+str(System_temp)+"K.txt",timestep,J2ev)
	#plot_heatflux(energy data2,Energy x axis min,Energy x axis max,i)
	TC.plot_heatflux(str(i)+"_Ener_equ_"+str(System_temp)+"K.txt",\
		Energy_xmin,Energy_xmax,i,Plot=True)

	#Calculate the temperature gradient in different way :dT/dL
	#temp_grad_dLdT(temperature profile data2,number of layers, number of fixed,number of bath,
	#temperature,temperature difference=20)
	TC.temp_grad_dTdL(str(i)+"_temp_equ_"+str(System_temp)+"K.txt",\
		number_layers,number_fixed,number_bath,System_temp)

	#计算热导率,默认使用拟合的温度梯度，若要使用dL/dT,需要将dLdT改成Ture
	TC.Thermal_conductivity("Thermal_conductivity.txt",thickness,i,dTdL=3)
'''
Thermal conductivities are saved in filename3   
dTdL=1,use fitting temperature gradient
dTdL=2,without including highest and lowest temperatures
dTdL=3,use directly temperature difference
'''
print('********************')
print('********************')
print('*****  Done!  ******')
print('********************')
print('********************')