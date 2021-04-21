import numpy as np 
import matplotlib.pyplot as plt 

def read_tc(tcfile):
	'''读取热导率文件，每行文件格式：x, y1, y2 y3
	返回x, y, y平均值, y的标准差'''
	tc = np.loadtxt(tcfile)
	# print(tc)
	x = tc[:,0]/32.06
	y = tc[:,1:4]
	y_mean = np.mean(y,axis=1)
	y_std = np.std(y,axis=1)
	# print(y_mean,y_std)
	return x, y, y_mean, y_std

def plt_tc(tc1,tc2=(1,1,1,1)):
	'''plot带误差棒热导率，'''

	x1, y1, y_mean1, y_std1 = tc1
	x2, y2, y_mean2, y_std2 = tc2

	plt.rc('font',family='Times New Roman',size=26)
	fig, ax = plt.subplots(figsize=(8,6))
	fig.subplots_adjust(bottom=0.2,left=0.2)

	# s1 = ax.errorbar(x1,y_mean1,yerr=y_std1,capsize=10,capthick=4,
	# 	fmt='ro:',mfc='white',mec='r',markersize=16,mew=2,alpha=0.4)
	# s2 = ax.errorbar(x2,y_mean2,yerr=y_std2,capsize=10,capthick=4,
	# 	fmt='b*:',mfc='white',mec='b',markersize=16,mew=2)

	s1 = ax.errorbar(x1,y_mean1,yerr=y_std1,capsize=10,capthick=4,
		fmt='ro:',mfc='white',mec='r',markersize=16,mew=2)
	s2 = ax.errorbar(x2,y_mean2,yerr=y_std2,capsize=10,capthick=4,
		fmt='b*:',mfc='white',mec='b',markersize=16,mew=2)
		
	ax.legend(handles=[s1,s2],labels=[
		r'$\mathregular{MoS_2}$/$\mathregular{MoS^{m}}_\mathregular{2}$',
		r'$\mathregular{MoS_2}$/$\mathregular{MoSe^{m}}_\mathregular{2}$']
		,loc='best', fontsize=22)

	ax.set_xlabel('Mass ratio (R)',fontsize=26,fontweight='bold')
	ax.set_ylabel('Thermal conductivity (W/m-K)',fontsize=26,fontweight='bold')
	ax.set_xticks([0,1,2,3,4,5,6])
	ax.set_yticks([0,5,10,15,20,25,30])
	plt.savefig('Thermal_conductivity.tiff',dpi=300,transparent=True)

	plt.show()



	return



if __name__ == '__main__':
	
	tcfile1 = './Thermal_conductivity_S.txt'
	tcfile2 = './Thermal_conductivity_Se.txt'
	tc1 = read_tc(tcfile1)
	tc2 = read_tc(tcfile2)
	'''若只有1条热导率点线图，
	只需：plt_tc(tc1),
	最多可画两条'''
	plt_tc(tc1,tc2)