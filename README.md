# photonic_sim
<img src="https://github.com/sopapado/photonic_sim/blob/master/photonic_sim%20logo.png" alt="drawing" width="200"/>

Photonic_sim is a collection of python methods for simulating propagation of fields in multilayered structures.
With photonic_sim you can calculate the electric and magnetic field of optical waves in arbitrary 1D photonic structures. The simulator assumes non-magnetic materials in current-free spaces.

## Description of method

The code is collecting all the boundary-condition related equations for electric and magnetic fields from every interface set by the user. Then it solves the system and manipulates the data in order to provide the field distribution along the structure for a given wavelength as well as transmission, reflection and loss spectrum for a given spectral region.

amp_distribution
field_distribution
calculate_spectrum
plot_distribution
plot_spetrum

There are some other methods that support the usage of these 5 methods but the user only needs the ones given in the list above. 

