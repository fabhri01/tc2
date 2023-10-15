#Importo librerias
import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from pytc2.sistemas_lineales import plot_plantilla

import warnings
warnings.filterwarnings('ignore')

# Para fijar el estilo de gráficos en el notebook
fig_sz_x = 10
fig_sz_y = 7
fig_dpi = 100
fig_font_size = 16

plt.rcParams['figure.figsize'] = (fig_sz_x, fig_sz_y)
plt.rcParams['figure.dpi'] = (fig_dpi)
plt.rcParams.update({'font.size':fig_font_size})

#Señal de ECG registrada a 1 kHz, con contaminación de diversos orígenes.
mat_struct = sio.loadmat('ecg.mat')

ecg_one_lead = mat_struct['ecg_lead']
ecg_one_lead = ecg_one_lead.flatten()
cant_muestras = len(ecg_one_lead)

fs = 1e3
nyq_frec = fs/2

#Plantilla
ws1 = 1.0 #Hz
wp1 = 3.0 #Hz
wp2 = 25.0 #Hz
ws2 = 35.0 #Hz

ripple = 0.5
atenuacion = 40

frecs = np.array([0.0,         ws1,         wp1,     wp2,     ws2,         nyq_frec   ]) / nyq_frec
gains = np.array([-atenuacion, -atenuacion, -ripple, -ripple, -atenuacion, -atenuacion])
gains = 10**(gains/20)

# -- FILTRO IIR --

iir = sig.iirdesign(wp=[wp1,wp2], ws=[ws1,ws2], gpass=ripple, gstop=atenuacion, analog=False, ftype='butter', output='sos', fs=1e3)

w, h = sig.sosfreqz(iir, worN = 2000, fs = fs)

#Modulo
plt.plot(w, 20*np.log10(abs(h)), label='IIR {:d} SOS'.format(iir.shape[0]))
plt.title('Filtro IIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Modulo [dB]')
plt.grid()
plt.axis([0, 100, -60, 5])
plot_plantilla(filter_type='bandpass', fpass=[wp1, wp2], ripple=ripple, fstop=[ws1, ws2], attenuation=atenuacion, fs=fs)
axes_hdl = plt.gca()
axes_hdl.legend()

#Fase
fase = np.angle(h)
plt.plot(w, fase)

plt.title('Fase filtro IIR')
plt.xlabel('Frecuencia [Hz]')
plt.ylabel('Fase [rad]')
plt.grid()
plt.axis([0, 100, -2*np.pi, 2*np.pi])

#Retardo de grupo
bn, an = sig.sos2tf(iir)
w, gd = sig.group_delay((bn, an), w=2000, whole=False, fs=fs)

plt.plot(w, gd)
plt.title('Retardo de grupo filtro IIR')
plt.xlabel('Retardo de grupo [muestras]')
plt.ylabel('Frecuencia [rad/muestras]')
plt.grid()

#Respuesta al impulso

imp = sig.unit_impulse(2000)
response = sig.sosfilt(iir, imp)

plt.plot(np.linspace(0, 1, 2000), response * 10)
plt.plot(np.linspace(0, 1, 2000), imp)

plt.title('Respuesta al impulso filtro IIR (x10)')
plt.xlabel('Tiempo [muestras]')
plt.ylabel('Amplitud')
plt.grid()

# Separar en 3 secciones distintas porque todo junto no va a graficar las 3 cosas

# Despues el resto es mas o menos igual a este codigo