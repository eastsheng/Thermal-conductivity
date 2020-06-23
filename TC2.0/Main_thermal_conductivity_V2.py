import thermal_conductivity_V2 as TCv2
##############################################
#case
k = 2
#系统尺寸(nm)
thickness = 0.33+0.35
#系统温度(K)
System_temp = 300
#层数
number_layers = 80
#固定层数
number_fixed = 1
#热浴层数
number_bath = 2
#单位换算
timestep=5e-7#ns
#热流方向 heatflux_direction=1：x方向为热流方向；heatflux_direction=2：y方向为热流方向
heatflux_direction = 2

#---------------------计算热导率从此处开始计算----------------------#
print(30*'-','Done!',30*'-')

for i in range(2,k+1):
	TCv2.Area('GRA_C3N_npt'+str(i)+'.data',"log.txt",i,thickness,heatflux_direction)

	#Read temperature profile for calculating temperature gradient
	#temp_grad(temperature profile data, number_layers,number_fixed,number_bath,i,fit_factor=3,Plot=False)
	TCv2.temp_grad(str(i)+"_temp_equ_"+str(System_temp)+"K.dat",\
		number_layers,number_fixed,number_bath,i,heatflux_direction,fit_factor=3,Plot=True)

	#Read input and output energies for calculating heat flux
	#heat_flux(energy data,timestep,J2ev)#J2ev=1.602763e-19#ev转换为J
	TCv2.heat_flux(str(i)+"_Ener_equ_"+str(System_temp)+"K.dat",i,timestep,Plot=True)
	'''   
	TempGrad_fator=1,use fitting temperature gradient.
	TempGrad_fator=2,without including highest and lowest temperatures,namely hot and cold bath.
	TempGrad_fator=3,use directly temperature difference.
	'''
	TCv2.Thermal_conductivity("Thermal_conductivity.txt",i,TempGrad_fator=1)

print(30*'-','Done!',30*'-')
