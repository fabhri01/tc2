#Importo librerias de la catedra
import sympy as sp
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2
from pytc2.general import print_latex

#Verificacion simbolica ej 2

s, z = sp.symbols('s z', complex = True)
k = sp.symbols('k', real = True, positive = True)

TB = 1/(s**2+s*np.sqrt(2)+1) #transferencia butterworth orden 2
FB = k*((z-1)/(z+1)) # transformada bilineal

Tz = sp.collect(sp.simplify(sp.expand(TB.subs(s, FB))), z)

print("Transferencia Butterworth orden 2 analógico")
print_latex('T(s) = ' + sp.latex(TB))

print("Transferencia del filtro digital correspondiente")
print_latex('T(z) = ' + sp.latex(Tz))

#como me pide repetir lo mismo para distintas frec defino una funcion
def graf(fc, fs):
    norma = fc #normalizo con fc
    fc = fc/norma
    fs = fs/norma

    #analogico
    z,p,k = sig.buttap(2)
    z,p,k = sig.lp2lp_zpk(z,p,k, wo=2*np.pi*fc)

    n, d = sig.zpk2tf(z,p,k)
    tf = sig.TransferFunction(n, d)
    
    #digital
    n_z, d_z = sig.bilinear(n, d, fs = fs)
    tf_z = sig.TransferFunction(n_z, d_z, dt = 1/fs)

    #printeo
    
    #polos y ceros
    tc2.pzmap(tf_z, annotations=True, filter_description="Butterworth digital de 2do orden", fig_id=2, digital=True)
    
    wrad_z, hh_z = sig.freqz(n_z, d_z)
    ww_z = wrad_z / np.pi
    
    plt.figure(1)

    #modulo
    plt.subplot(2,1,1)
    plt.grid(visible=True)
    plt.title("Respuesta de módulo")
    plt.ylabel("Módulo [dB]")
    plt.plot(ww_z, 20*np.log10(abs(hh_z)), color = 'b', label='Butterworth digital de 2do orden')

    axes_hdl = plt.gca()
    axes_hdl.legend()

    #fase
    plt.subplot(2,1,2)
    plt.grid(visible=True)
    plt.title("Respuesta de fase")
    plt.ylabel("Fase [deg]")
    plt.xlabel("Frecuencia normalizada [fs/2]")
    plt.plot(ww_z, np.angle(hh_z, deg=True), color = 'b', label='Butterworth digital de 2do orden')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    
graf(1e3, 100e3)



