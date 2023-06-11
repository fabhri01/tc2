import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
#libreria de la catedra
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys, pretty_print_SOS, tf2sos_analog
from pytc2.general import print_subtitle

#Filtro pasabajos prototipo
#wp = 1
ws = 4
alpha_max = 1
alpha_min = 30

eps = np.sqrt(10**(alpha_max/10) - 1)
orden = np.ceil(0.5*(np.log10((10**(alpha_min*0.1) - 1)/ eps**2) / np.log10(ws)))
z, p, k = sig.buttap(orden)
num, den = sig.zpk2tf(z,p,k)
print_subtitle("Transferencia normalizada filtro pasabajos Butterworth prototipo")
sos_lp = tf2sos_analog(num, den)
sos_lp[sos_lp < 1e-6] = 0.0
pretty_print_SOS(sos_lp)
analyze_sys(sos_lp)

#Filtro pasa altos
num, den = sig.lp2lp(num, den, eps**(-1/orden)) #primero desnormalizo el filtro pasa bajos por usar wbutter
num, den = sig.lp2hp(num, den) #transformo el filtro pasabajos de MP desnormalizado a un filtro pasa altos
sos_hp = tf2sos_analog(num, den)
print_subtitle("Transferencia filtro pasa altos MP")
pretty_print_SOS(sos_hp)
analyze_sys(sos_hp)
