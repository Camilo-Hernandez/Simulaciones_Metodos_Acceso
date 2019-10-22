## Importar librerias
from scipy.stats import poisson, binom, expon	# Funciones que guardan todas las propiedades de la VA Poisson, Exponencial y Binomial
from math import exp
import matplotlib.pyplot as plt 	# Graficador

## ------------ ALOHA ranurado ------------ ##
'''
Cuando fueron generadas tramas nuevas en el slot anterior, los nodos deciden transmitir,
con un 100% de probabilidad, en el slot actual
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
'S'=P(1) es la eficiencia entendida como la fraccion de los slots que contienen una Tx exitosa: S = g*exp(-g)
'''
slots = 5				# Slots totales a analizar
print("Slots: ",slots)
## Eficiencia S en funcion de la probabilidad de retransmision qr
# Parametros
lamb = [x/10 for x in range(7,12)]				# Vector de tasas medias de generacion de tramas
tramas = [list(poisson.rvs(mu=i,size=slots)) for i in lamb]     # Tramas generadas aleatoriamente por slot por cada tasa del vector lamb
qr = [x/100 for x in range(1,101,25)]                           # Vector de probabilidades de retransmision de los nodos BL
n = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Cantidad inicial de nodos BL por slot para cada qr y cada lambda
iBL = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Porciones iniciales de los n que intentan retransmitir en cada slot
g = [[[0 for x in range(slots)] for y in range(len(qr))] for z in range(len(lamb))] # Tasa de Tx (funcion de n)
Exitos = [[0 for x in range(len(qr))] for y in range(len(lamb))] # Slots o tramas con Tx exitosa
S = [[0 for x in range(len(qr))] for y in range(len(lamb))] # Eficiencia usando tasa de Tx (funcion de qr y lambda)
Efc = [[0 for x in range(len(qr))] for y in range(len(lamb))] # Eficiencia usando porcentaje de tramas transmitidas (funcion de qr y lambda)

# Ciclo que calcula las eficiencias
for k in range(len(lamb)):
    tram=tramas[k]
    print("Media Poisson de tramas/slot: ",lamb[k])
    print("Tramas/slot: ", tram)
    for i in range(len(qr)):
        for j in range(1,slots):
            iBL[k][i][j] = binom.rvs(n=n[k][i][j-1],p=qr[i])
            if tram[j-1]+iBL[k][i][j]>1:
                n[k][i][j]=n[k][i][j-1]+tram[j-1]                 # n aumenta en la cantidad de tramas generadas en el anterior debido a colision
            elif tram[j-1]==1 and iBL[k][i][j]==0:
                Exitos[k][i] += 1
                n[k][i][j]=n[k][i][j-1]                           # +1 slot con Tx exitosa. Y n no cambia
            elif tram[j-1]==0 and iBL[k][i][j]==1:
                Exitos[k][i] += 1
                n[k][i][j]=n[k][i][j-1]-1                         # +1 slot con Txexitosa. Y n disminuye en 1
            elif tram[j-1]==0 and iBL[k][i][j]==0:
                n[k][i][j]=n[k][i][j-1]

#Efc = [[100*y*sum(k) for y,k in zip(x,tramas)] for x in Exitos]
#g = [[lam + n * m for n in x] for x,m in zip(n,qr)]     # g = constante lam + matriz n .* vector qr
#S = [100*sum(x*exp(-x) for x in row)/slots for row in g]

plt.plot(qr,S)
plt.plot(qr,Efc)
plt.title("Eficiencias segun la probabilidad de retransmision de los nodos backlogged para un lambda dado")
plt.xlabel("Probabilidad qr")
plt.ylabel("S(g(n)) - Efc(tramas enviadas)")
plt.legend(["S= g(n)*exp{-g(n)}","Efc = Enviadas/Generadas"])
plt.show()

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