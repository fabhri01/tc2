from scipy.signal import TransferFunction
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as ss

from pytc2.sistemas_lineales import analyze_sys, pretty_print_bicuad_omegayq, tf2sos_analog, pretty_print_SOS

from pytc2.general import print_subtitle

#----------------------Definicion de las variables---------------------

xi = 0.349
n = 3
alpha_max = 0.5

#---------------------------Definicion de H-------------------------

[z, p, k] = ss.cheb1ap(n, alpha_max)
[num, den] = ss.zpk2tf(z, p, k)
num,den = ss.lp2bp(num,den, wo=1, bw=0.2)

sos_bp = tf2sos_analog(num, den)
sos_bp[sos_bp < 1e-6] = 0.0

#-----------------------------Visualizacion----------------------------

print_subtitle('Transferencia factorizada del Filtro pasa banda Cheby de tercer orden')
pretty_print_SOS(sos_bp)
print_subtitle('Transferencia factorizada y parametrizada del Filtro pasa banda Cheby de tercer orden')
pretty_print_SOS(sos_bp, mode='omegayq')

analyze_sys(sos_bp, "Cheby")

  
