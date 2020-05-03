import numpy as np 
import matplotlib.pyplot as plt 
import pprint
from tabulate import tabulate
from tqdm import tqdm
import w2rgb as w2rgb
import addcopyfighandler as copyfig
import my_functions.my_functions as mf
from scipy.signal import find_peaks as find_peaks 
from scipy.interpolate import interp1d as interp1d
import sys
import mpmath as mp
from IPython.lib.pretty import pretty
from mpmath import fp
import time
import scipy.constants as sc
import cmath

def transform_ns(lamda,ns):
    
    # this pretty function is checking whether there is spectally defined refractive index and takes the value for the current wavelength and puts it in the initial ns
    # it transforms the ns in a simple 1D array
    # IF ns is not to be transformed meaning that the ns is simple (single values of indices) 
    # it returns the ns. 

    ns_new = ns[:]
    for i in range(len(ns)):
        if np.size(ns[i])>1:
            fmat = interp1d( ns[i][0], ns[i][1] )
            ns_new[i] = complex(fmat(lamda))
            
    return ns_new


def field_distribution(lamda, ns, ds, thetai, y0, amps, dx,pol): 

    # lamda : the wavelength in m that you want the distribution to be calculated
    # ns : the refractive indices for every layer (including input and output region) 
    # ds : the thicknesses in every layer (len(ds) == len(ns))
    # amps : an array with all the aplitudes of roward and backward propagating waves in every layer (amp_distribution function is giving this)
    # dx : the resolution of the plot in m
    # returns field, nx, x_plot. 
    # field: is the complex field distribution along the structure with dx resolution 
    # nx : is the refractive index distribution given by the user with dx resolution this time ( to be plotted alond with the field)
    # x_plot : is the x-axis value for plotting field vs distance x. 
    
    # Transform ns 
    ns_new = transform_ns(lamda, ns)

    
    xtemp = 0
    xs = []

    for d in ds:
        xtemp = d + xtemp
        xs.append(xtemp)
    
    k = 2*np.pi/lamda
    N = len(ds)-1

    xs = list(xs)
    xs2 = [0]+xs
    
    
    Ex = []
    Ey = []
    Ez = []
    Hx = []
    Hy = []
    Hz = []

    x_plot = []
    nx = []
    y = y0

    # a sweet loop that runs in every layer and adds together forward and backward waves for estimating the distribution inside it.
    
    for i in range(N+1):
        
        exp = np.exp
        n = np.int(np.abs(ds[i])/dx)
        x_layer = np.linspace(xs2[i],xs2[i+1],n)
        #x_layer = mp.matrix(x_layer)
        nx = nx + list(np.ones((n,1))*ns_new[i])
        x_plot = x_plot + list(x_layer)
        
        #this commented lines are not used - I leave it here to remember how to do that
        #c0 = np.asarray(list(map(mp.exp,1j*ns_new[i]*k*x_layer)))
        #c1 = np.asarray(list(map(mp.exp,-1j*ns_new[i]*k*x_layer)))

        k1 = 1j*(x_layer*np.cos(thetai[i])+y*np.sin(thetai[i]))*k*ns_new[i]
        k2 = 1j*(-x_layer*np.cos(thetai[i])+y*np.sin(thetai[i]))*k*ns_new[i]

        if pol == 'p':

            Ec0x = -np.sin(thetai[i])*exp(k1)
            Ec1x = -np.sin(thetai[i])*exp(k2)
            Ec0y = np.cos(thetai[i])*exp(k1)
            Ec1y = -np.cos(thetai[i])*exp(k2)
            Ec0z = np.zeros(n)
            Ec1z = np.zeros(n)

            Hc0x = np.zeros(n)
            Hc1x = np.zeros(n)
            Hc0y = np.zeros(n)
            Hc1y = np.zeros(n)
            Hc0z = -ns_new[i]*exp(k1)
            Hc1z = -ns_new[i]*exp(k2)

            layer_efieldx = amps[2*i]*Ec0x + amps[2*i+1]*Ec1x
            layer_efieldy = amps[2*i]*Ec0y + amps[2*i+1]*Ec1y
            layer_efieldz = amps[2*i]*Ec0z + amps[2*i+1]*Ec1z

            layer_hfieldx = amps[2*i]*Hc0x + amps[2*i+1]*Hc1x
            layer_hfieldy = amps[2*i]*Hc0y + amps[2*i+1]*Hc1y
            layer_hfieldz = amps[2*i]*Hc0z + amps[2*i+1]*Hc1z
                 
            lefx = list(layer_efieldx)
            lefy = list(layer_efieldy)
            lefz = list(layer_efieldz)

            lhfx = list(layer_hfieldx)
            lhfy = list(layer_hfieldy)
            lhfz = list(layer_hfieldz)

            Ex = Ex + lefx 
            Ey = Ey + lefy
            Ez = Ez + lefz 

            Hx = Hx + lhfx 
            Hy = Hy + lhfy
            Hz = Hz + lhfz 
        
        if pol == 's' : 

            Ec0x = np.zeros(n)
            Ec1x = np.zeros(n)
            Ec0y = np.zeros(n)
            Ec1y = np.zeros(n)
            Ec0z = exp(k1)
            Ec1z = exp(k2)

            Hc0x = -ns_new[i]*exp(k1)*np.sin(thetai[i])
            Hc1x = -ns_new[i]*exp(k2)*np.sin(thetai[i])
            Hc0y = ns_new[i]*exp(k1)*np.cos(thetai[i])
            Hc1y = -ns_new[i]*exp(k2)*np.cos(thetai[i])
            Hc0z = np.zeros(n)
            Hc1z = np.zeros(n)
            
            layer_efieldx = amps[2*i]*Ec0x + amps[2*i+1]*Ec1x
            layer_efieldy = amps[2*i]*Ec0y + amps[2*i+1]*Ec1y
            layer_efieldz = amps[2*i]*Ec0z + amps[2*i+1]*Ec1z

            layer_hfieldx = amps[2*i]*Hc0x + amps[2*i+1]*Hc1x
            layer_hfieldy = amps[2*i]*Hc0y + amps[2*i+1]*Hc1y
            layer_hfieldz = amps[2*i]*Hc0z + amps[2*i+1]*Hc1z
                 
            lefx = list(layer_efieldx)
            lefy = list(layer_efieldy)
            lefz = list(layer_efieldz)

            lhfx = list(layer_hfieldx)
            lhfy = list(layer_hfieldy)
            lhfz = list(layer_hfieldz)

            Ex = Ex + lefx 
            Ey = Ey + lefy
            Ez = Ez + lefz 

            Hx = Hx + lhfx 
            Hy = Hy + lhfy
            Hz = Hz + lhfz 
            

    Ex = list(map(complex, Ex))
    Ey = list(map(complex, Ey))
    Ez = list(map(complex, Ez))
    Hx = list(map(complex, Hx))
    Hy = list(map(complex, Hy))
    Hz = list(map(complex, Hz))

    x_plot = list(map(float, x_plot)) 
    nx = list(map(complex, nx))


    

    return np.asarray(Ex), np.asarray(Ey), np.asarray(Ez), np.asarray(Hx), np.asarray(Hy), np.asarray(Hz), np.asarray(nx), np.asarray(x_plot)


def plot_distribution(E,H, nx, x, lamda, filename ):
    
    rgb = w2rgb.wavelength_to_rgb(lamda*1e9, gamma = 1)
    plt.rcParams.update({'font.size': 18})
    fig,ax1 = plt.subplots(figsize=(10.0,4.0))
    ax1.plot(np.asarray(x)/1e-9,E, color = rgb, linewidth = 3,label = '$|E|^2$')
    ax1.plot(np.asarray(x)/1e-9,H, color = 'gray', linewidth = 2,linestyle = '--',label = '$|H|^2$')
    ax1.set_ylim(0)
    plt.legend()
    ax2 = ax1.twinx()
    nx = np.squeeze(nx)
    ax2.fill_between(np.asarray(x)/1e-9,np.squeeze(np.abs(nx)),np.zeros(len(nx))+np.min(np.abs(nx))-0.1, color = 'cy')
    ax2.set_ylim(np.min(np.abs(nx))-0.1,np.max(np.abs(nx)+0.5))
    ax1.set_xlabel('distance [nm]')
    ax1.set_ylabel('normalized [a.u.]')
    ax2.set_ylabel('refractive index', color = (0,0.8,0.8))
    plt.xlim(x[0]/1e-9, x[-1]/1e-9)
    plt.title('$\lambda : $'+str(np.around(lamda*1e9,3))+ ' nm')
    plt.tight_layout()
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.patch.set_visible(False)
    
    plt.savefig(filename+'.png')
    mf.save_fig(fig, filename)
    plt.clf()

 
def calculate_spectrum(lamda, ns, ds, l0, l1, dlamda, pol, theta_init, linit, amp_init, x0, y0):



    # lamda : wavelength in which the amplitudes are calculated
    # ns : 1D array gives consecutive values for all the layers including input and output media. At least size 3.
    # ds : 1D array of consecutive thickness values of the layers. Same size as ns. At least size 3.
    # spectrum_span is the window of spectrum to be calculated around the lamda set.
    # dlamda is the resolution in m that you want your spectrum to be computed
    # returns R, T and lamdas
    # R is the intensity reflectivity value 
    # T is the intensity transmission value
    # lamdas is a 1D array with the wavelengths values

    Rspec = []
    Tspec = []
 
    
    # for i in range(len(ns)):
    #     if np.size(ns[i])>1:
    #         lamdas = np.linspace( ns[i][0][0], ns[i][0][-1], np.int( (ns[i][0][-1]-ns[i][0][0])/dlamda ) )
    
    # lamdas = np.linspace(lamda-spectrum_span/2,lamda+spectrum_span/2,np.int(spectrum_span/dlamda))
    
    lamdas = define_lamdas( lamda, dlamda, l0,l1, ns )
    
    xtemp = x0
    xs = []
    for d in ds:
        xtemp = d + xtemp
        xs.append(xtemp)
    
    k = 2*np.pi/lamda
    N = len(ns)-2
    xs = list(xs)
    xs2 = [0]+xs

    print('# CALCULATING SPECTRUM #')
    ns = [ns[0]] + ns + [ns[-1]]
    ds = [0.1e-9] + ds + [0.1e-9]
    exp = np.exp
    for lamda in tqdm(lamdas):
        
        
        ns_new = transform_ns(lamda, ns)

        
        amps, xs, thetai = amp_distribution(lamda, ns, ds, pol, theta_init, linit, amp_init, x0, y0)
        # the reflection and transmission of the structure will be given by the amplitudes of in first and last layer
        R = np.abs(amps[1])**2
        T = np.abs(amps[-2])**2*np.abs((ns_new[-1]/ns_new[0]))
        

        Rspec.append(R)
        Tspec.append(T)

    print('')
    Tspec = list(map(float,Tspec))
    Rspec = list(map(float,Rspec))

    return np.asarray(Rspec), np.asarray(Tspec), lamdas


def define_lamdas(lamda, dlamda, l0,l1, ns ):
    
    # this method is using the ns array where the user has put any dispersion data
    # of the material (complex refractive index VS wavelength) and defines the 
    # array that is going to be used for the calculation of the spectrum.
    # IF no material dispersion is given meaning that ns is a simple 1D array with
    # single elements then the spectrum_span and lamda given by user will be used
    # to define the array lamdas for the spectrum calculation.

    ind = []
    simple = True


    for i in range(len(ns)):
        if np.size(ns[i])>1: # if the element of ns is not size 1 then it assuses dispersion data
            ind.append(i)
            simple = False

    l_min = []
    l_max = []
    dl = [] 

    for i in ind:
        l_min.append(ns[i][0][0])
        l_max.append(ns[i][0][-1]) 
        dl.append(ns[i][0][-1] - ns[i][0][0])

    l_min = [l0] + l_min
    l_max = [l1] + l_max
    

    lamda_min = np.max(l_min)
    lamda_max = np.min(l_max)


    lamdas = np.linspace(lamda_min, lamda_max, np.int((lamda_max-lamda_min)/dlamda))



    return lamdas 

def amp_distribution(lamda, ns, ds,pol, theta_init, linit, amp_init, x0, y0):
            
    k = 2*np.pi/lamda
    N = len(ds)-1 # number of interfaces

    #checking that the inputs ns and ds are same 
    # size. Return error if not
    if len(ns)!=len(ds) : 
        print(" ERROR : ns and ds not same size")
        exit()

    A = np.zeros((2*N,2*N+2)) + 0j
    exp= np.exp

    # Transform ns
    ns_new = transform_ns(lamda, ns)

    # using the thickness value given from user to define the x-coordinates of the interfaces (xs)
    xtemp = x0
    xs = []
    for d in ds:
        xtemp = d + xtemp
        xs.append(xtemp)


    # calculate the angles of propagation for all layers

    thetai = np.zeros(N+1) + 0j
    thetai[linit] = theta_init

    for i in range(N-linit):
        nl = ns_new[i+linit] # refractive index on the left of interface
        nr = ns_new[i+1+linit] # refractive index on the right of interface
        thetai[i+linit+1] =   cmath.asin(nl*np.sin(thetai[i+linit])/nr)

    for i in range(linit):
        nl = ns_new[linit-i] # refractive index on the left of interface
        nr = ns_new[linit-i-1] # refractive index on the right of interface
        thetai[linit-i-1] =  cmath.asin(nl*np.sin(thetai[i+linit])/nr)


    for i in range(N):
        x = xs[i] # coordinate of interface
        nl = ns_new[i] # refractive index on the left of interface
        nr = ns_new[i+1] # refractive index on the right of interface
        thetal = thetai[i]
        thetar = thetai[i+1]
        y = y0

    # field equations that enter in matrix A (check notes)

        k1l = 1j*(x*np.cos(thetal)+y*np.sin(thetal))*k*nl
        k2l = 1j*(-x*np.cos(thetal)+y*np.sin(thetal))*k*nl
        k1r = 1j*(x*np.cos(thetar)+y*np.sin(thetar))*k*nr
        k2r = 1j*(-x*np.cos(thetar)+y*np.sin(thetar))*k*nr

        if pol == 'p' :

            a11 = np.cos(thetal)*exp(k1l)
            a12 = -np.cos(thetal)*exp(k2l)
            a13 = -np.cos(thetar)*exp(k1r)
            a14 = np.cos(thetar)*exp(k2r)

            a21 = nl*exp(k1l)
            a22 = nl*exp(k2l)
            a23 = -nr*exp(k1r)
            a24 = -nr*exp(k2r)

            A[2*i,2*i] = a11 
            A[2*i,2*i+1] = a12  
            A[2*i,2*i+2] = a13
            A[2*i,2*i+3] = a14

            A[2*i+1,2*i] = a21 
            A[2*i+1,2*i+1] = a22
            A[2*i+1,2*i+2] = a23
            A[2*i+1,2*i+3] = a24

        if pol == 's' :

            a11 = exp(k1l)
            a12 = exp(k2l)
            a13 = -exp(k1r)
            a14 = -exp(k2r)

            a21 = np.cos(thetal)*nl*exp(k1l)
            a22 = -np.cos(thetal)*nl*exp(k2l)
            a23 = -np.cos(thetar)*nr*exp(k1r)
            a24 = np.cos(thetar)*nr*exp(k2r)

            A[2*i,2*i] = a11 
            A[2*i,2*i+1] = a12  
            A[2*i,2*i+2] = a13
            A[2*i,2*i+3] = a14

            A[2*i+1,2*i] = a21 
            A[2*i+1,2*i+1] = a22
            A[2*i+1,2*i+2] = a23
            A[2*i+1,2*i+3] = a24    

    # now we have a y = Ax  system to solve 
    # but we need to take away two columns of initiation
    # we form then a square matrix Asq


    Asq = A[:,1:-1]
    if linit != 0 : s = -(A[:,2*linit] + A[:,2*linit+1] )
    if linit == 0 : s = -A[:,0]
    sols = np.linalg.solve(Asq,s)

    if linit == 0 : amps = [amp_init] + list(sols) + [0]
    if linit != 0 : 
        amps = [0] + list(sols) + [0]
        amps[2*linit] = amps[2*linit]+amp_init/2
        amps[2*linit+1] = amps[2*linit+1]+amp_init/2

    return amps, xs, thetai

def plot_spectrum(T, R, lamdas, filename ):
    
   
    # T = np.asarray(T)
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.0,4.0))
    plt.grid(zorder = 20)
    plt.plot(lamdas/1e-9, T,'darkred',linewidth = 2)
    plt.fill_between(lamdas/1e-9, T,np.zeros(len(T)),color = 'darkred',linewidth = 2, zorder = 10, label = 'T')
    plt.plot(lamdas/1e-9, R+T ,'k',linestyle = '--', label = 'T + R')
    plt.xlim(lamdas[0]/1e-9,lamdas[-1]/1e-9)
    plt.xlabel('wavelength [nm]')
    plt.ylabel('normalized intenstiy')
    #plt.legend(loc = 'best')#,  bbox_to_anchor=(1.5, 1, 0, 0) )
    plt.ylim(0,1+0.05)
    plt.tight_layout()
    
    plt.savefig(filename+'T.png')
    mf.save_fig(fig, filename+ 'T')
    plt.clf()
    
    plt.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots(figsize=(6.0,4.0))
    plt.grid(zorder = 20)
    plt.plot(lamdas/1e-9, R,'darkgreen',linewidth = 2)
    plt.fill_between(lamdas/1e-9, R,np.zeros(len(T)),color = 'darkgreen',linewidth = 2, zorder = 10, label = 'R')
    # plt.plot(lamdas/1e-9, R+T ,'k',linestyle = '--', label = 'T + R')
    plt.xlim(lamdas[0]/1e-9,lamdas[-1]/1e-9)
    plt.xlabel('wavelength [nm]')
    plt.ylabel('normalized intensity')
    # plt.legend(loc = 'best')#,  bbox_to_anchor=(1.5, 1, 0, 0) )
    plt.ylim(0,1+0.05)
    plt.tight_layout()
    
    plt.savefig(filename+'R.png')
    mf.save_fig(fig, filename+ 'R')
    plt.clf()


def print_peak_wavelengths(T, lamdas):
    
    # this function gives the wavelength values of the peaks in transmission spectrum

    peaks, prop = find_peaks(T)
    
    lpeak = [list(np.linspace(1,len(peaks),len(peaks)))]+[list(lamdas[peaks]/1e-9)]
    print('### PEAK WAVELENGTHS [nm] ###' )
    print(tabulate(np.transpose(lpeak)))
    return peaks