

# se importan librerías y funciones de la cátedra

from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot
##%matplotlib qt5


#Funciones de aproximacion por maxima planicidad
def MP_aprox(ws, alpha_max, alpha_min):

    xi = ( np.sqrt( 10**(alpha_max / 10) - 1 ))    #Por lo general w_p=1    

    order=np.ceil( 0.5 * ( np.log10( ( 10**(alpha_min * 0.1) - 1 ) / xi**2)  / np.log10(ws) ) ) 

    z,p,k = sig.buttap (order)
    N, D = sig.zpk2tf(z, p, k)
    N, D = sig.lp2lp(N, D, xi**(-1/order) ) #Esta función cambia el filtro para otra frecuencia. Es la renormalización del filtro para epsilon distinto de 1
    
    z,p,k=sig.tf2zpk(N,D)
    
    return N, D


#Simulación de la función transferencia diseñada

#wp=1 // siempre se normalizará a 1
ws=2
alfa_max=1
alfa_min=12

N,D = MP_aprox(ws, alfa_max, alfa_min)
T = TransferFunction(N, D)

bodePlot(T, fig_id=1)
pzmap(T, fig_id=2 )
GroupDelay(T, fig_id=3)