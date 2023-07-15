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

