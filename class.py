import numpy as np
import matplotlib.pyplot as plt

def EulerRichardson(T, dt,N, R,I, A):
	'''
	Funcion para solucionar ecuaciones diferenciales utilizando
	el metodo de Euler-Richardson
	T  : tiempo de integracion 			 [int]
	dt : step del metodo       	     	 [float]
	N  : numero de personas    			 [int]
	R  : numero de recuperados iniciales [int]
	I  : numero de infectados iniciales  [int] 
	A  : proporcion de infectados versus casos totales confirmados [float]
	'''



	# Constantes
	wk = 7
	yr = 365
	Tc = 5       # tiempo promedio en lka que la gente se junta
	Tr = 14      # Tiempo necesario para que alguien se recupere
	M = 0.001    # Mortalidad
	delta = 0.05
	alpha = A
	gamma = 1/Tr
	beta = 1/Tc

	T = T*wk

	# Arrays para guardar los datos 
	General = []       # Lista de todos los datos
	Rf = np.zeros(T)   # Recovered people
	If = np.zeros(T)   # Infected people
	Sf = np.zeros(T)   # Suceptible people
	x  = np.zeros(T)   # Counter 
	Npp = np.zeros(T)  # S + I + R
	R2 = np.zeros(T)   # People death
	R1 = np.zeros(T)   # People recovered
	C = np.zeros(T)    # Beds needed
	P = np.zeros(T)    # Positive Test taken
	
	# Condiciones iniciales
	Rf[0] = R
	If[0] = I
	Sf[0] = N-If[0]
	Npp[0] = N
	R1[0] = 0
	R2[0] = 0
	C[0] = 0

	# Loop para el metodo Euler-Richardson

	for n in range(0,T-1):
	    Sfmid = Sf[n] + 0.5*dt*(-beta*(If[n]*Sf[n])/N)
	    Rfmid = Rf[n] + 0.5*dt*gamma*If[n]
	    Ifmid = If[n] + 0.5*dt*((beta*If[n]*Sf[n])/N - gamma*If[n])
	    
	    Sf[n+1] = Sf[n] + dt*(-beta*(Ifmid*Sfmid)/N)
	    Rf[n+1] = Rf[n] + dt*gamma*Ifmid
	    If[n+1] = If[n] + dt*((beta*Ifmid*Sfmid)/N - gamma*Ifmid)
	    
	    
	    Npp[n+1] = Sf[n+1] + Rf[n+1] + If[n+1]
	    R2[n+1] = M*Rf[n]
	    R1[n+1] = Rf[n] - R2[n]
	    P[n+1] = If[n]/alpha + P[n]
	    C[n+1] = delta*If[n]

	    x[n+1] = n

	General.append(Sf)
	General.append(Rf)
	General.append(If)
	General.append(R1)
	General.append(R2)
	General.append(Npp)
	General.append(P)
	General.append(C)
	General.append(x)


	return General

def Plots(Sf,Rf,If,R1,R2,Npp,P,C,x):

	plt.figure(figsize = (15,6))
	plt.plot(x,Rf, label = 'Recovered')
	plt.plot(x,Sf, label = 'Suceptible')
	plt.plot(x,If, label = 'Infected')
	plt.plot(x,Npp, label = 'S+I+R', ls = '--')
	plt.legend(loc = 'best', fontsize = 13)
	plt.xlabel('Days', labelpad = 15, fontsize = 15)
	plt.ylabel('Pople number ', labelpad = 15, fontsize = 15)
	plt.title('S-I-R epidemic model', fontsize = 15)
	#plt.savefig('./SIR_Eu_Rich.png')
	plt.show()


	plt.figure(figsize = (20,5))
	plt.subplot(2,2,1)
	plt.plot(x,C, label= 'Beds ', color = 'crimson')
	plt.legend(loc = 'best', fontsize = 15)
	plt.xlabel('Days', labelpad = 15, fontsize = 15)
	plt.ylabel('Beds number ', labelpad = 15, fontsize = 15)

	plt.subplot(2,2,2)
	plt.plot(x,P, label= 'Positive test (R1)', color = 'k')
	plt.legend(loc = 'best', fontsize = 15)
	plt.xlabel('Days', labelpad = 15, fontsize = 15)
	plt.ylabel('Possitive tests', labelpad = 15, fontsize = 15)


	plt.subplot(2,2,3)
	plt.plot(x,R2, label= 'People Death (R2)', color = 'magenta')
	plt.legend(loc = 'best', fontsize = 15)
	plt.xlabel('Days', labelpad = 15, fontsize = 15)
	plt.ylabel('Pople number ', labelpad = 15, fontsize = 15)


	plt.subplot(2,2,4)
	plt.plot(x,R1, label= 'People recovered (R1)', color = 'olive')
	plt.legend(loc = 'best', fontsize = 15)
	plt.xlabel('Days', labelpad = 15, fontsize = 15)

	#plt.savefig('./beds_positive_test_200I.png')
	plt.show()


Data = EulerRichardson(156, 0.11,100000, 0,1, 3)

Plots(Data[0],Data[1], Data[2], Data[3],Data[4], Data[5], Data[6],Data[7], Data[8])


