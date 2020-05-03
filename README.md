# photonic_sim
<img src="https://github.com/sopapado/photonic_sim/blob/master/photonic_sim%20logo.png" alt="drawing" width="200"/>

Photonic_sim is a collection of python methods for simulating propagation of fields in multilayered structures.
With photonic_sim you can calculate the electric and magnetic field of optical waves in arbitrary 1D photonic structures. The simulator assumes non-magnetic materials in current-free spaces.

## Description of method

The code is collecting all the boundary-condition related equations for electric and magnetic fields from every interface set by the user. Then it solves the system and manipulates the data in order to provide the field distribution along the structure for a given wavelength as well as transmission, reflection and loss spectrum for a given spectral region.

## Getting started 

Let's start with an example. Assuming that we want to simulate a multilayer structure consisting of three layers with different refractive indices sitting in air as shown in the image (n: refractive index). We assume also that the optical field is p-polarized and has an incident angle of 0 degrees.



To set the geometric and optical parameters we use the following code:

```python

import photonic_sim.photonic_sim as ps

lamda = 500e-9
theta = 0
pol = 'p'
ns = [1, 1.5, 1.2, 1.5, 1]
ds = [200e-9, 125e-9, 250e-9, 125e-9, 200e-9]

```
- **lamda** : the wavelength for which we will calculate the distribution.
- **theta** : the angle of incidence of the optical field at the structure.
- **pol** : polarization mode of the incoming field. Possible string values 's' and 'p'.
- **ns** : the values of the refractive indices including the air at the beggining and at the end of the structure.
- **ds** : the thicknesses of the layers. NOTE that the first and last values can be arbitrary because they refer to the air but they should always be included.

As you see the air parts of the structure are treated as additional layers in our geometry. This means that you can assume that you enter or exit the structure with any kind of material you want. Note that there is no reflection at the end and at the start of the structure. 
We are all set for our first calculation. We will use the *amp_distribution* function. 

```python

amps, xs, thetai = ps.amp_distribution(lamda, ns, ds,pol, theta, linit = 0, amp_init = 1, x0 = 0, y0 = 0):

```
(for now I do not comment on linit, amp_init, x0 and y0 variables )

This function solves our mathematical problem and gives us:
- **amps**: a 1D array with the amplitudes of the fields propagating right (A) and propagating left (B) in every layer i (i = 0:number_of_layers-1). For example A1 is the right propagating wave amplitude in layer 1 and B3 is the left propagating wave in layer 3. Given that notation the amplitude list will have the form amps : [A0, B0, A1, B1, A2, B2, ... ]
- **xs**: the x-coordinates of every interface (computed simply by the *ds* array given by the user)
- **thetai** : the angles of propagation in every layer



