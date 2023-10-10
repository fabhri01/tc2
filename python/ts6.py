#Importo librerias de la catedra
import sympy as sp
import pytc2.cuadripolos as tc2
from pytc2.general import print_latex

#Verificacion simbolica ej 2

s = sp.symbols('s', complex = True)
Z1, Z3 = sp.symbols('L1 L3', complex = True)
Y1 = 1/(s*Z1)
Y2 = s*sp.symbols('C2', complex = True)
Y3 = 1/(s*Z3)
R = sp.symbols('R', real = True, positive = True)
G = 1/R

#MAI 

MAI = sp.Matrix([
                    [Y1, -Y1, 0, 0],
                    [-Y1, Y1+Y2+Y3, -Y2, -Y3],
                    [0, -Y2, Y2 + G, -G],
                    [0, -Y3, -G, Y3 + G]
                ])

con_detalles = True

#Verifico transferencia
tf = tc2.calc_MAI_vtransf_ij_mn(MAI, 2, 3, 0, 2, verbose = con_detalles)

print('Cálculo de la transferencia')
print_latex('T(s) = ' + sp.latex(tf))

#Reemplazo valores para obtenerla normalizada
tf_n = sp.simplify(tf.subs(Z1, 3/2))
tf_n = sp.simplify(tf_n.subs(Y2, s*(4/3)))
tf_n = sp.simplify(tf_n.subs(Z3, 1/2))
tf_n = sp.simplify(tf_n.subs(R, 1))

print_latex('T(s) = ' + sp.latex(tf_n))

#Verifico impedancia de entrada
con_detalles = False
Ze = tc2.calc_MAI_impedance_ij(MAI, 0, 2, verbose = con_detalles)
print('Cálculo de la impedancia de entrada')
print_latex('Ze(s) = ' + sp.latex(Ze))