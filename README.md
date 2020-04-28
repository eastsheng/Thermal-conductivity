#***************************Thermal-conductivity***************************#                                                               
Calculating thermal conductivity;                                                                                                         
Just to urge myself to insist on writing code.                                                                                             

1. The function that size(data,number_layers) is used to read the system size for calculating thermal conductivity. Wherein 'numer_layer'
is the number of layers of system by dividing in NEMD simulation.

2. The function that temp_grad(tempdata1,writetempdata2,number_layers) is used to read and write the temperature profile file for plotting and fitting, obtaining the temperature gradient finally.

3. The function that temp_grad_dLdT(tempdata1,number_layers,number_fixed,number_bath) is a another way to calculate the temperature gradient, that is, dL/dT---x (heat flux) direction size of system/temperature difference without including fixed and bath layers.

4. heat_flux(energydata1,writeenergydata2,timestep,J2ev). As the name suggests, this is a function that read energy of input and output. Then, the energies is dealed to obtain a avage energy profile for obtaining the heat flux.

5. plot_temp(writetempdata2,layer_fixed,number_fixed,number_bath,i). Plotting and fitting the temperature profile, then, a temperature gradient is obtained. 'i' is the cases number.

6. plot_heatflux(filename2,Energy_xmin,Energy_xmax,i). Plotting and fitting. Energy_xmin and Energy_xmax are min and max value about time by fitting, respectively.

7. A function to calculate thermal conductivity, Thermal_conductivity(ThermalConductivityResultFile, system thickness,i,dLdT=False), if you want to use '3' to calculate temperature gradient, you should modify to dLdT=True.
