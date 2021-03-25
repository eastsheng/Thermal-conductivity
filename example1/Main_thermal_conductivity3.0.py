import ThermalConductivity as TC

##############################################
#case
k = 1
#系统尺寸(nm)
thickness = 0.61
#系统温度(K)
System_temp = 300
#层数
number_layers = 120
#固定层数
number_fixed = 2
#热浴层数
number_bath = 2
#单位换算
timestep=5e-7#ns
#热流方向 heatflux_direction=1：x方向为热流方向；heatflux_direction=2：y方向为热流方向
heatflux_direction = 1

path_tc = '../example1/'

#---------------------计算热导率从此处开始计算----------------------#
print(30*'-','Done!',30*'-')
tc = TC.ThermalConductivity()
for i in range(1,k+1):
	relax_data = path_tc+'MoS2_NPT.data'#'./0811/TC3.0/MoS2_NPT.data'
	tc.read_size(relax_data,i,thickness,heatflux_direction)

	#Read temperature profile for calculating temperature gradient
	temperaturefile = str(i)+"_temp_equ_"+str(System_temp)+"K.dat"#'./0811/TC3.0/'+str(i)+"_temp_equ_"+str(System_temp)+"K.dat"
	tc.temp_grad(path_tc,temperaturefile,number_layers,number_fixed,number_bath,fit_factor=3,Plot=True)

	#Read input and output energies for calculating heat flux

	heatfluxfile =str(i)+"_Ener_equ_"+str(System_temp)+"K.dat" #'./0811/TC3.0/'+str(i)+"_Ener_equ_"+str(System_temp)+"K.dat"
	tc.heat_flux(path_tc,heatfluxfile,timestep)

	'''   
	TempGrad_fator=1,use fitting temperature gradient.
	TempGrad_fator=2,without including highest and lowest temperatures,namely hot and cold bath.
	TempGrad_fator=3,use directly temperature difference.
	'''
	result = 'Thermal_conductivity.txt'#"./0811/TC3.0/Thermal_conductivity.txt"
	tc.thermal_conductivity(result,TempGrad_fator=1)

	logname = 'log.txt'#'./0811/TC3.0/log.txt'
	tc.logfile(logname)

print(30*'-','Done!',30*'-')
