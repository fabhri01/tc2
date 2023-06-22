# importo librerias
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
#libreria de la catedra
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys, pretty_print_SOS, tf2sos_analog
from pytc2.general import print_subtitle

#Enunciado y datos del gráfico de módulo y fase
fc = 300
fz = 100

wc = fc/fc
wz = fz/fc

orden = 3

#Pasabajos prototipo
Omega_c = 1/wc
Omega_z = 1/wz

n_z = [0, 1/(wz**2)] #cero de transmisión que agregare al array num

z, p, k = sig.buttap(orden)
num, den = sig.zpk2tf(z, p, k) #obtengo num y den de un butter orden 3

n_z = np.append(num, n_z) #asi tengo un numerador s²+wz² para el cero de transmision

n_z, den = sig.lp2lp(n_z * (1/Omega_z**2), den) #para 0dB

sos_lp = tf2sos_analog(n_z, den) 
sos_lp[sos_lp < 1e-6] = 0.0
print_subtitle("Transferencia prototipo pasabajos Butterworth orden 3")
pretty_print_SOS(sos_lp, mode='omegayq')

num2, den2 = sig.lp2hp(n_z, den)
sos_hp = tf2sos_analog(num2, den2)
sos_hp[sos_hp < 1e-6] = 0.0
print_subtitle("Transferencia objetivo pasa altos Butterworth")
pretty_print_SOS(sos_hp, mode='omegayq')

analyze_sys(sos_hp)