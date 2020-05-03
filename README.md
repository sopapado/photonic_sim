# photonic_sim
<img src="https://github.com/sopapado/photonic_sim/blob/master/photonic_sim%20logo.png" alt="drawing" width="200"/>

Photonic_sim is a collection of python methods for simulating propagation of fields in multilayered structures.
With photonic_sim you can calculate the electric and magnetic field of optical waves in arbitrary 1D photonic structures. The simulator assumes non-magnetic materials in current-free spaces.

## Description of method

The code is collecting all the boundary-condition related equations for electric and magnetic fields from every interface set by the user. Then it solves the system and manipulates the data in order to provide the field distribution along the structure for a given wavelength as well as transmission, reflection and loss spectrum for a given spectral region.

## Getting started 

Let's start with an example. Assuming that we want to simulate a multilayer structure consisting of three layers with different refractive indices sitting in air as shown in the image (n: refractive index)


To set the geometric and optical parameters we use the following code:

```markdown

import photonic_sim.photonic_sim as ps

lamda = 500e-9
theta = 0
ns = [1, 1.5, 1.2, 1.5, 1]
ds = [200e-9, 125e-9, 250e-9, 125e-9, 200e-9]

```
- **lamda** : the wavelength for which we will calculate the distribution.
- **theta** : the angle of incidence of the optical field at the structure.
- **ns** : the values of the refractive indices including the air at the beggining and at the end of the structure.
- **ds** : the thicknesses of the layers. NOTE that the first and last values can be arbitrary because they refer to the air but they should always be included.



amp_distribution
field_distribution
calculate_spectrum
plot_distribution
plot_spetrum

There are some other methods that support the usage of these 5 methods but the user only needs the ones given in the list above. 

