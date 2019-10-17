## Importar librerias
from scipy.stats import poisson, expon, binom	# Funciones que guardan todas las propiedades de la VA Poisson, Exponencial y Binomial
from random import random, expovariate
from math import exp
from numpy import cumsum		# Creador de vectores acumulativos
import matplotlib.pyplot as plt 	# Graficador

##  Aloha ranurado  ##
'''
Cuando se generan tramas nuevas, los nodos deciden transmitir con un 100% de probabilidad en el siguiente slot
Si en el sgte slot de ser generadas no se transmiten o hay colision,
se convierten en backlogged los nodos y decidiran transmitir con un 19% de probabilidad en el siguiente slot
Se le llama intento BL cuando un BL decide retransmitir
Los slots tienen duracion de 1s
'lamb' es una lista que guarda la cantidad de tramas generadas por slot (VA Poisson con tasa media mu=0.2 tramas por segundo por ejemplo)
	Ej: Si 'lamb' igual a [0, 0, 0, 0, 1, 0,] indica que se analizaron 6 slots y en el quinto slot se genero una trama, en el resto no
'lambBL' es una variable entera que guarda la cantidad de tramas no enviadas en el sgte slot de ser generadas (tramas de nodos BL). Es acumulativa.
'lambBL' es lo que se llamaria 'n' en el mapa de Markov modelado por la cantidad de nodos BL
'lambBL' crece en 1 o mas cuando hay colision, que sucede en un slot cuando
	1. Se generan mas de una trama
	2. Hay mas de un intento BL
	3. Se genera una trama y un intento BL (o mas de cada uno)
'lambBL' no cambia en el sgte slot si:
	1. No hay tramas generadas ni intentos BL
	2. No se generan tramas nuevas pero hay colision entre los nodos BL
	3. Se envia una trama nueva exitosamente sin que hayan intentos BL
'lambBL' decrece en 1 cuando no se generan tramas nuevas y hay un unico intento BL
'qr' es una lista de VA Binomiales que indican la cantidad de intentos BL por slot donde: n=lambBL, p=probabilidad de retransmitir
'g' es la tasa de transmision por slot = lamb + qr
	Donde 'g' == 1, en ese slot hubo una transmision exitosa, y si 'qr' == 1 en ese mismo slot, lambBL -= 1
	Donde 'g' > 1, en ese slot hubo una colision y lambBL += tramas generadas en el anterior slot
	Donde 'g' < 1 que pasa????
'S' es la eficiencia total del sistema, es una lista de porcentajes por slot. S = g*exp(-g)
lambBL PODRIA SER UNA LISTA ACUMULATIVA POR SLOTS
'''
nodos = 10				# Nodos totales
slots = 10				# Slots totales
lambBL = 0				# Cantidad de nodos BL inicial
muqr = 0.2				# Tasa media de retransmision
mulamb = 1				# Tasa media de generacion de tramas
lamb = list(poisson.rvs(mu=mulamb,size=slots))
qr=[0]*slots
for i in range(len(lamb)):
	qr[i] += binom.rvs(p=muqr,n=lambBL)	# Intento de retransmision en el sgte slot
	if lamb[i] == 0 or lamb[i] == 1: colision = 0;
	elif lamb[i] > 1:
		colision = 1
		lambBL += lamb[i]	# Aumento de nodos BL

qr[1]=binom.rvs(n=lambBL,p=muqr)
qr = list(binom.rvs(n=12, p=0.1, size=slots))
g = list(map(lambda x,y: x+y, lamb,qr))
qr = bernoulli.rvs(p=0.2)		# Los nodos backlogged en promedio deciden transmitir en el 20% de los slots siguientes
P0 = exp(-g)*100			# Probabilidad de no tener intentos de transmision en un slot
P1 = g*exp(-g)*100			# Probabilidad de tener un unico intento de transmision o de retransmision
Pc = 100-P0-P1				# Probabilidad de colision
S.append(P1)				# Lista de probabilidades de transmitir
nlist.append(n)				# Lista de nodos backlogged
glist.append(g)				# Lista de tasas de transmision
slots.append(slots[-1]+1)
rango = cumsum([P0,P1,Pc])
coin = random()*rango[2]
if P1 == rango[rango>coin][0]: n -= 1	# Disminucion de 1 nodo backlogged al darse una transmision exitosa
n += lamb				# Crecimiento indefinido de n nodos backlogged
slots.pop()

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
n = 0
t = expon.rvs(1, size=n)		# Tiempo de espera para intentar retransmitir
x = 1/t					# Tasa de retransmision de un nodo
lamb = poisson.rvs(mu=1)		# Tasa de generacion de tramas nuevas en la red
G = lamb + x				# Tasa de intentos de transmision en el sistema en el estado n (Poisson process)
Pe = exp(-2*G)				# Probabilidad de exito
n += lamb				# Incremento de n nodos backlogged
h=1