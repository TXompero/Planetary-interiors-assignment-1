import numpy as np
import matplotlib.pyplot as plt
from burnman import Composite, minerals, Layer, PerplexMaterial, Planet, BoundaryLayerPerturbation
from burnman.tools.chemistry import formula_mass

# Planet layers definition
R_planet = 2440*1.e3
alpha = 0.6895
beta = 0.984

# Determine inner core material and parameters

core_material = minerals.SE_2015.hcp_iron()

core_radius = 300
core = Layer('inner core',radii = np.linspace(0,R_planet*alpha,100))
core.set_material(core_material)
core.set_temperature_mode('adiabatic',
                          temperature_top = 2000)

core.set_pressure_mode(pressure_mode = 'self-consistent',
                       pressure_top = 11.65e9,
                       gravity_bottom = 0.1)

core.make()

mantle_material = minerals.SLB_2024.olivine([0.333,0.667])
mantle = Layer('mantle',radii = np.linspace(R_planet*alpha,R_planet*beta,100))
mantle.set_material(mantle_material)

tbl_perturbation = BoundaryLayerPerturbation(radius_bottom=R_planet*alpha,
                                             radius_top=R_planet*beta,
                                             rayleigh_number=1.e4,
                                             temperature_change=1000.,
                                             boundary_layer_ratio=0.5)

#mantle.set_temperature_mode('perturbed-adiabatic',
#                            temperatures = tbl_perturbation.temperature(np.linspace(R_planet*alpha,R_planet*beta,100)))

mantle.set_temperature_mode('user-defined',
                            temperatures = np.linspace(2000,800,100))

crust_material = minerals.SLB_2024.quartz()
crust = Layer('crust',radii = np.linspace(R_planet*beta,R_planet,10))
crust.set_material(crust_material)
crust.set_temperature_mode('user-defined',
                           np.linspace(800,440,10))

mercury = Planet('Mercury',[core, mantle, crust],verbose=True)
mercury.make()




fig = plt.figure(figsize=(7, 3), dpi = 300)
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]
ax[0].plot(mercury.density,mercury.radii/1.e3)
ax[1].plot(mercury.temperature,mercury.radii/1.e3)
for i in range(1):
    ax[i].set_ylim(mercury.radii[0]/1.e3,
                   mercury.radii[-1]/1.e3)
    ax[i].set_ylabel('Radius (km)')

ax[0].set_xlabel('Density (kg/m$^3$)')
ax[1].set_xlabel('Temperature (K)')


fig.set_layout_engine('tight')
fig.savefig("figures/burnman_results.png")

print(f'mass = {mercury.average_density:.3e}')
print(f'moment of inertia factor= {mercury.moment_of_inertia_factor:.4f}')