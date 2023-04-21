## Comparison of results

- [LAMMPS website result](https://docs.lammps.org/compute_heat_flux.html): 
  - 0.29 W/mK
- Running LAMMPS locally automatically yields the following results: 
  - average conductivity: 0.382412490910382[W/mK] @ 70 K, 0.0257443666020476 /A\^3
- The result obtained by using script `greenkubo_tc_from_HFACF.py` integral is:
  - The average thermal conductivity is calculated as 0.3894952057937662 (W/(mK)), error is 0.0
- Using the script `greenkubo_tc_from_HeatFlux.py` autocorrelation and integration gets the result:
  - kappa = 0.2909348733644516 (W/mK), error = ± 0.00024125366861261992