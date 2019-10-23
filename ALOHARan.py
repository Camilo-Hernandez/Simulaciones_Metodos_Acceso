'''
Este programa ealiza graficas de la eficiencia del ALOHA ranurado
en funcion de la probabilidad de retransmision de los nodos backlogged
y la tasa media de generacion de tramas en la red
2 eficiencias se calculan: S y Efc
'''
## Importar librerias
from scipy.stats import poisson, binom	# Funciones que guardan todas las propiedades de la VA Poisson y Binomial
from math import exp
import pylab as pl; import numpy as np	# Para graficas en 3D y manejo de vectores y matrices
from mpl_toolkits.mplot3d import Axes3D  # Para hacer graficas en el espacio 3D.


## ------------ ALOHA ranurado (descripcion del programa) ------------ ##
'''
Cuando fueron generadas tramas nuevas en el slot anterior, los nodos deciden transmitir,
con un 100% de probabilidad, en el slot actual.
Si en el slot anterior hubo colision,
significa que aparecieron nuevos nodos backlogged igual a la cantidad de tramas generadas en el anterior,
de los cuales algunos decidiran transmitir en el slot actual con probabilidad 'qr'
Los slots tienen duracion de 1 unidad de tiempo. Estan representados por las posiciones de las sub-listas mas internas
formadas en las iteraciones del ciclo 'for' mas interno o anidado
'tramas' es VA Poisson con tasa media 'lamb', y es lista bi-dimensional que indica
las tramas generadas en cada slot por su respectiva 'lamb'
	Ej: Si 'tramas' igual a [0, 0, 0, 1, 0,] indica que se analizaron 5 slots y en el 4to slot se genero una trama, en el resto no
'n' es la variable que modela el mapa de Markov del ALOHA ranurado
'n' es una lista tri-dimensional que guarda la cantidad total de nodos BL por slot, para cada 'qr' en el rango
y para cada trama generada aleatoriamente por su respectiva tasa 'lamb' diferente (elementos por slot, vectores por 'qr' y matrices por 'lamb')
Se le llama intento BL cuando un nodo BL decide retransmitir
'iBL' es una lista de VA Binomiales tri-dimensional que indican la porcion de 'n' que intenten retransmitir en el slot actual,
para cada 'qr' y 'lamb' en sus rangos. Los parametros de la distribucion son los intentos n='n' y la probabilidad p='qr'
'n' crece en el slot actual respecto del anterior en 1 o mas cuando hay colision, es decir si:
	1. Se generan mas de una trama en el slot anterior
	2. Hay mas de un intento BL en el slot actual
	3. La suma del numero de tramas generadas en el slot anterior e intentos BL en el actual es mayor que 1 (engloba los dos anteriores)
'n' no cambia de valor en el slot actual respecto del anterior si:
	1. No hay tramas generadas en el anterior ni intentos BL en el actual
	2. No se generaron tramas nuevas en el anterior pero si hay mas de 1 intento BL en el actual
	3. Se genera un unica trama nueva en el slot anterior y no hay intentos BL en el actual (Tx exitosa)
'n' decrece en 1 si no se generaron tramas nuevas en el anterior y hay un unico intento BL en el actual (Tx exitosa)
'g' es la tasa de intentos de transmision de tramas por slot: g(n) = lamb + n*qr
'g' es una lista tri-dimensional: contiene matrices por cada 'lamb', vectores por cada 'qr', tasas por cada 'slot'
'S'=P(1) es la eficiencia total del sistema entendida como la fraccion de los slots que contienen una Tx exitosa.
'S' es una lista bi-dimensional (o matriz) que guarda los promedios de los vectores que hay por cada 'qr' de la lista 'g'. S = promedio(g*exp(-g) por cada 'qr')
'Efc' es la eficiencia total del sistema entendida como la razon entre de tramas enviadas y las generadas.
'Efc' y 'S' son bi-dimensionales: las filas (o vectores/sub-listas) son por cada 'lamb', y las columnas (o elementos de cada sub-lista) son por cada 'qr'
El programa crea las graficas en el espacio 3D de ambas eficiencias en funcion de los vectores 'qr' y 'lamb' que varian entre 0.1-1 y 0.1 y 2, respectivamente
'''

## ------------ Algoritmo ------------ ##

# Parametros
slots = 5							# Slots totales a analizar
print("Slots: ",slots)
lamb = [x/10 for x in range(1,21)]				# Vector de tasas medias de generacion de tramas
tramas = [list(poisson.rvs(mu=i,size=slots)) for i in lamb]     # Tramas generadas aleatoriamente por slot por cada tasa del vector lamb
qr = [x/100 for x in range(1,101,2)]                           	# Vector de probabilidades de retransmision de los nodos BL
n = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Cantidad nodos BL por slot para cada qr y cada lambda
iBL = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Porciones iniciales de los n que intentan retransmitir por slot
g = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Tasa de Tx. Es funcion de lambda, n y qr.
S = [[0 for x in range(len(qr))] for y in range(len(lamb))]  	 # Eficiencia usando promedio(g). Es funcion de lambda y qr
Exitos = [[0 for x in range(len(qr))] for y in range(len(lamb))] # Slots o tramas con Tx exitosa
Efc = [[0 for x in range(len(qr))] for y in range(len(lamb))]    # Eficiencia usando porcentaje de tramas transmitidas respecto a las generadas. Es funcion de lambda y qr

# Ciclo que recorre la lista 'tramas'
for k in range(len(lamb)):
    tram=tramas[k]
    print("Media Poisson de tramas/slot: ",lamb[k])
    print("Tramas/slot: ", tram)
    # Ciclo que recorre la lista 'qr'
    for i in range(len(qr)):
	# Ciclo que calcula 'n' y 'Exitos' por cada slot
        for j in range(1,slots):
            iBL[k][i][j] = binom.rvs(n=n[k][i][j-1],p=qr[i])	# Nodos BL que quieren retransmitir
            if tram[j-1]+iBL[k][i][j]>1:			# Condicion 3 para que 'n' crezca en 1 o mas
                n[k][i][j]=n[k][i][j-1]+tram[j-1]		# (tambien es la 2 para que 'n' no cambie)
            elif tram[j-1]==1 and iBL[k][i][j]==0:		# Condicion 3 para que 'n' no cambie
                Exitos[k][i] += 1				# +1 slot con Tx exitosa
                n[k][i][j]=n[k][i][j-1]                        
            elif tram[j-1]==0 and iBL[k][i][j]==1:
                Exitos[k][i] += 1				# +1 slot con Txexitosa. Y n disminuye en 1
                n[k][i][j]=n[k][i][j-1]-1                       
            elif tram[j-1]==0 and iBL[k][i][j]==0:		# Condicion 1 para que 'n' no cambie
                n[k][i][j]=n[k][i][j-1]

Efc = [[y1/sum(x2) for y1 in x1] for x1,x2 in zip(Exitos,tramas)] # Tramas enviadas/tramas generadas
g = [[[lam + q * nBL for nBL in nx] for nx,q in zip(nlam,qr)] for lam, nlam in zip(lamb,n)]
S = [[100*sum(x*exp(-x) for x in row)/slots for row in glmatrix] for glmatrix in g] # S = promedio(g*exp(-g))

## ------------ Graficas de las eficiencias en funcion de los vectores 'lamb' y 'qr' ------------ ##

# Superficie S
fig = pl.figure(); ax = Axes3D(fig)
X,Y=np.meshgrid(np.array(lamb),np.array(qr))
Z = np.array(S).T
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=pl.cm.winter)
ax.set_title("Eficiencia S del ALOHA Ranurado")
ax.set_xlabel('Tasa media lambda', fontsize=15)
ax.set_ylabel('Probabilidad qr',fontsize=15)
ax.set_zlabel('S = %prom[g(n)*exp{-g(n)}]',fontsize=11)
pl.show()

# Superficie Efc
fig = pl.figure(); ax = Axes3D(fig)
Z = np.array(Efc).T
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=pl.cm.summer)
ax.set_title("Eficiencia Efc")
ax.set_xlabel('Tasa media lambda', fontsize=15)
ax.set_ylabel('Probabilidad qr',fontsize=15)
ax.set_zlabel('Efc = %Enviadas/Generadas',fontsize=15)
pl.show()
