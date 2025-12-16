# import ThermalConductivity as TC
import readlammpsdata as rld
import ReadLammpsTraj as rlt
import fastdataing.fastdataing as fd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

def cal_area(lmp,hfd="z"):
	"""
	Calculating the cross sectional area
	lmp: relax lammps data
	hfd: heatflux direction
	
	"""
	lx = rld.read_len(lmp,direction="x")
	ly = rld.read_len(lmp,direction="y")
	lz = rld.read_len(lmp,direction="z")
	area_map = {
		"x": ly * lz,
		"y": lx * lz,
		"z": lx * ly,
	}
	area = area_map[hfd]*1e-20 # A^2 --> m^2
	print(f"cross sectional area: {area} m^2")
	return lx,ly,lz,area


def filter_temp_data(temp, threshold=20):
    data = temp[:, [0, 3]]
    col2 = data[:, 1]
    vmax, vmin, vmean = col2.max(), col2.min(), col2.mean()
    mask = (col2 <= vmax + threshold) & (col2 >= vmax - 20 - threshold)
    return data[mask]

def temp_grad(tempfile,ff=3):
	"""
	Calculating the temperature gradient 
	tempfile: temperature distribution file
	ff: fit_factor
	"""
	temp = rlt.average_avechunk(tempfile)
	temp[:, 0] *= (lz * 0.1) / temp.shape[0]
	filtered_data = filter_temp_data(temp)
	z,temp = filtered_data[:,0],filtered_data[:,1]
	fx,fy,f_fit = fd.polyfitting(z[ff:-ff],temp[ff:-ff],degree=1)
	tg = f_fit[1]
	print(f"temp_grad: {tg} K/nm")
	with plt.style.context(['science', 'no-latex']):
		fig, ax = plt.subplots()
		ax.scatter(z, temp, c="b")
		ax.plot(fx,fy, c="r")
	ax.autoscale(tight=True)
	ax.set(
	    xlabel=r"$z$ (nm)",
	    ylabel=r"$T$",
	    xlim=(-1, 11),
	    ylim=(185, 215)
		)
	fig.tight_layout()
	# fig.savefig(f"cn_N7_Zn2+.png", dpi=300)
	plt.show()
	plt.close()
	return tg


def heatflux(energyfile,ff=3):
	"""
	Calculating the heat flux density
	energyfile: energy file
	ff: fit_factor
	"""
	kcalpermol2J = 4186.0/6.02214e23

	energy = np.loadtxt(energyfile,skiprows=1).reshape(-1,3)
	t = energy[:,0]*timestep*1e-6 # ns
	E1,E2 = energy[:,1]*kcalpermol2J,energy[:,2]*kcalpermol2J
	q = (abs(E1)+abs(E2))/2.0
	qx,qy,Q_fit = fd.polyfitting(t,q,degree=1)
	hf = abs(Q_fit[1]) # J/ns
	print(f"heat flux: {hf} J/ns")
	with plt.style.context(['science', 'no-latex']):
		fig, ax = plt.subplots()
		ax.scatter(t, E1*1e17, c="b")
		ax.scatter(t, E2*1e17, c="g")
		# ax.plot(fx,fy, c="r")
	ax.autoscale(tight=True)
	ax.set(
	    xlabel=r"$t$ (ns)",
	    ylabel=r"Energy ($\times 10^{-17} (J)$)",
	    # xlim=(-1, 11),
	    ylim=(-5, 5)
		)
	fig.tight_layout()
	# fig.savefig(f"cn_N7_Zn2+.png", dpi=300)
	plt.show()
	plt.close()

	return hf



if __name__ == '__main__':
	k = 3
	#系统温度(K)
	T = 200
	timestep=1 # fs
	# heatflux_direction
	hfd = "z"
	path = './200K/'
	ff = 5

	#--------------------- Start ----------------------#
	# tc = TC.ThermalConductivity()
	for i in range(1,k+1):
		relax_data = path+f'{k}_nvt_relax_{T}.data'
		# print(help(rld))
		lx,ly,lz,area = cal_area(lmp=relax_data,hfd=hfd)
		tempfile = f"{path}{i}_temp_equ_{T}K.dat"
		energyfile = f"{path}{i}_Ener_equ_{T}K.dat"
		tg = temp_grad(tempfile,ff=ff)
		hf = heatflux(energyfile,ff=ff)

		kappa = -hf/(tg*area)
		print(f"kappa = {kappa} W/mK")
