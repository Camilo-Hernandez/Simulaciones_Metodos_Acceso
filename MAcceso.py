## Importar librerias
from scipy.stats import poisson, binom, expon	# Funciones que guardan todas las propiedades de la VA Poisson, Exponencial y Binomial
from random import random, expovariate
from math import exp
from numpy import cumsum		# Creador de vectores acumulativos
import matplotlib.pyplot as plt 	# Graficador

## ------------ ALOHA ranurado ------------ ##
'''
Cuando fueron generadas tramas nuevas en el slot anterior, los nodos deciden transmitir con un 100% de probabilidad en el slot actual
Si en el slot anterior hubo colision,
significa que aparecieron nuevos nodos backlogged igual a la cantidad de tramas generadas en el anterior,
de los cuales algunos decidiran transmitir en el slot actual con probabilidad 'qr'
'n' es la variable que modela el mapa de Markov del ALOHA ranurado
'n' es una lista que guarda la cantidad total de nodos BL por slot, o mejor dicho, de # tramas que no han sido enviadas luego del slot en que fueron generadas
Los slots tienen duracion de 1 unidad de tiempo. Estan representados por las iteraciones del ciclo 'for'
'tramas' es la lista que indica la cantidad de tramas generadas en cada slot (VA Poisson con tasa media lamb)
	Ej: Si 'tramas' igual a [0, 0, 0, 1, 0,] indica que se analizaron 5 slots y en el 4to slot se genero una trama, en el resto no
Se le llama intento BL cuando un nodo BL decide retransmitir
'iBL' es una lista de VA Binomiales que indican la porcion de 'n' que intenten retransmitir en el slot actual,
donde: n='n', p='qr' (no confundir la variable 'n' con n que es uno de los dos parametros de la VA Binomial)
'n' crece en 1 o mas cuando hay colision, que sucede en el slot actual cuando
	1. Se generan mas de una trama en el slot anterior
	2. Hay mas de un intento BL en el slot actual
	3. La suma del numero de tramas generadas en el slot anterior e intentos BL en el actual es mayor que 1 (engloba los dos anteriores)
'n' no cambia de valor en el slot actual si:
	1. No hay tramas generadas en el anterior ni intentos BL en el actual
	2. No se generaron tramas nuevas en el anterior pero si hay mas de 1 intento BL en el actual
	3. Se genera un unica trama nueva en el slot anterior y no hay intentos BL en el actual
'n' decrece en 1 si no se generaron tramas nuevas en el anterior y hay un unico intento BL en el actual
'g' es la tasa de intentos de transmision de tramas por slot: g = lamb + n*qr
'Efc' es la eficiencia total del sistema entendida como la razon entre de tramas enviadas y las generadas
'S' es la eficiencia entendida como la fraccion de los slots que contienen una Tx exitosa: S = g*exp(-g)
'''
slots = 15				# Slots totales a analizar
Exitos = 0                              # Slots o tramas con Tx exitosa
n = [0]*slots			        # Cantidad inicial de nodos BL por slot (lista de ceros)
qr = 0.2                                # Probabilidad de retransmision (constante)
iBL = [0]*slots                         # Porciones iniciales de los n que intentan retransmitir
lamb = 0.5				# Tasa media de generacion de tramas (constante)
g = [0]*slots                           # Tasa inicial de Tx (varia por slot y es funcion de n)
tramas = list(poisson.rvs(mu=lamb,size=slots))  # (En Python, la media se denota 'mu' en la funcion 'poisson.rvs()')
# Ciclo que calcula 'n' e 'iBL' en el slot actual
# 'iBL' es aleatorio segun los 'n' que hayan acumulados hasta el slot anterior
# 'n' se calcula segun las 'tramas' el anterior slot y los 'iBL' en el slot actual
for i in range(slots):
    iBL[i] = binom.rvs(n=n[i-1],p=qr)
    if tramas[i-1]+iBL[i]>1:
        n[i]=n[i-1]+tramas[i-1]                 # n aumenta en la cantidad de tramas generadas en el anterior debido a colision
    elif (tramas[i-1]==1 and iBL[i]==0):
        Exitos += 1                              # +1 slot con Tx de trama exitosa. Y n no cambia
    elif tramas[i-1]==0 and iBL[i]==1:
        Exitos += 1
        n[i]=n[i-1]-1	                        # +1 slot con Tx de trama exitosa. Y n disminuye en 1

Efc = Exitos / sum(tramas)
g = [x+lamb for x in binom.mean(n,qr)]
P1 = [x*exp(-x)*100 for x in g]		        # Probabilidad de tener un unico envio de trama (exito)
print('De ',sum(tramas),' tramas generadas, se enviaron ',Exitos,'. La eficiencia del sistema fue del ',Efc,'%.', sep='')

## Grafica de barras
#width = 0.6  # the width of the bars

#fig, ax = plt.subplots()
##rects1 = ax.bar(slots - width/2, S, width, label='Men')
#rects2 = ax.bar(slots + width/2, nlist, width, label='Women')

## Add some text for labels, title and custom x-axis tick labels, etc.
##ax.set_ylabel('Scores')
#ax.set_title('Aloha ranurado con 30 nodos')
##ax.set_xticks(x)
#ax.set_xticklabels(slots)
#ax.legend()


## Aloha no ranurado
'''Si ocurre una colision,los nodos involucrados
intentaran transmitir la trama despues de un tiempo aleatorio'''
#n = 0
#t = expon.rvs(1, size=n)		# Tiempo de espera para intentar retransmitir
#x = 1/t					# Tasa de retransmision de un nodo
#tramas = poisson.rvs(mu=1)		# Tasa de generacion de tramas nuevas en la red
#G = tramas + x				# Tasa de intentos de transmision en el sistema en el estado n (Poisson process)
#Pe = exp(-2*G)				# Probabilidad de exito
#n += tramas				# Incremento de n nodos backlogged
#h=1