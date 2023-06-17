# importo librerias
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
#libreria de la catedra
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys, pretty_print_SOS, tf2sos_analog
from pytc2.general import print_subtitle

pi = np.pi

#D = 80e-6
#norma_w = 1/D
##solo hace falta el orden

orden = 3

z, p, k = sig.besselap(orden, norm='delay')

num, den = sig.zpk2tf(z, p, k)
sos_lp = tf2sos_analog(num, den)
sos_lp[sos_lp < 1e-6] = 0.0
pretty_print_SOS(sos_lp, mode='omegayq')
analyze_sys(sos_lp)

#Cuentas sintetización circuital
#componentes del integrado y ecuaciones de diseño
R1 = 50e3
C1 = 1000e-12

C2 = C1
R2 = R1 
R4 = R1

#Rf1 = 1 
#Rf2 = 2

#w_n = np.sqrt(R2/(R1*C1*C2*Rf1*Rf2))
#Q = ( (1 + (R4*(Rg+Rq)) / (Rg*Rq)) / (1 + (R2/R1)) ) * np.sqrt( (R2*C1*Rf1)/(R1*C2*Rf2) )
#k = (1 + (R1/R2)) / (Rg * (1/Rg + 1/Rq + 1/R4))

#Desarrollo y calculo de componentes externos

w_o = 2.542 #lo que dio el numpy
Q = 0.691

f = 2e3
norma_w = 2*pi*f
w_n = (w_o * norma_w)
print(f'frecuencia desnormalizada: {w_n}')

Rg = 50e3
Rq = 30.8e3 #variable
