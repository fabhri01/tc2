# importo librerias
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
#libreria de la catedra
from pytc2.sistemas_lineales import pzmap, GroupDelay, bodePlot, analyze_sys, pretty_print_SOS, tf2sos_analog
from pytc2.general import print_subtitle

n = 3
alpha_max = 0.5
BW = 1/5
w_o_n = 1

#%matplotlib qt5

z, p, k = sig.cheb1ap(n, alpha_max)
num, den = sig.zpk2tf(z, p, k)
sos_lp = tf2sos_analog(num, den)
sos_lp[sos_lp < 1e-6] = 0.0
print_subtitle("Transferencia prototipo pasabajos Chebyshev orden 3")
pretty_print_SOS(sos_lp, mode='omegayq')

num2, den2 = sig.lp2bp(num, den, w_o_n, BW)
sos_bp = tf2sos_analog(num2, den2)
sos_bp[sos_bp < 1e-6] = 0.0
print_subtitle("Transferencia objetivo pasa altos Chebyshev")
pretty_print_SOS(sos_bp, mode='omegayq')

analyze_sys(sos_bp)