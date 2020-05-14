import numpy as np
import matplotlib.pyplot as plt

class EulerRichardson():
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
	def __init__(self, T, dt ,N, R, I, A):
		"""
			Inicializa constantes, listas y condiciones iniciales.

			TO DO: Acuerdate de ser consistente con las variables. Tienes que añadirle el self si no
			no van a existir en este metodo.
		"""

		wk = 7
		yr = 365

		self.T = T*wk
		self.dt = dt
		self.N = N
		self.R = R
		self.I = I
		self.A = A


		Tc = 5
		Tr = 14
		self.M = 0.001
		self.delta = 0.05
		self.alpha = self.A
		self.gamma = 1/Tr
		self.beta = 1/Tc


		self.rf = np.zeros(self.T)   
		self.if_ = np.zeros(self.T)   
		self.sf = np.zeros(self.T)   
		self.x  = np.zeros(self.T)   
		self.npp = np.zeros(self.T)  
		self.r2 = np.zeros(self.T)   
		self.r1 = np.zeros(self.T)   
		self.c = np.zeros(self.T)    
		self.p = np.zeros(self.T)    
	
		
		
		# Condiciones iniciales
		self.rf[0] = self.R
		self.if_[0] = self.I
		self.sf[0] = self.N-self.if_[0]
		self.npp[0] = self.N
		self.r1[0] = 0
		self.r2[0] = 0
		self.c[0] = 0

	def solve(self):
		"""
			Loop para el metodo Euler-Richardson

			TO DO: Acuerdate de ser consistente con las variables. Tienes que añadirle el self si no
			no van a existir en este metodo.
		"""
		for n in range(0,self.T-1):
		    sfmid = self.sf[n] + 0.5*self.dt*(-self.beta*(self.if_[n]*self.sf[n])/self.N)
		    rfmid = self.rf[n] + 0.5*self.dt*self.gamma*self.if_[n]
		    if_mid = self.if_[n] + 0.5*self.dt*((self.beta*self.if_[n]*self.sf[n])/self.N - self.gamma*self.if_[n])
		    
		    self.sf[n+1] = self.sf[n] + self.dt*(-self.beta*(if_mid*sfmid)/self.N)
		    self.rf[n+1] = self.rf[n] + self.dt*self.gamma*if_mid
		    self.if_[n+1] = self.if_[n] + self.dt*((self.beta*if_mid*sfmid)/self.N - self.gamma*if_mid)
			#if_[n+1] = if_[n] + self.dt*((beta*if_mid*sfmid)/self.N - gamma*if_mid)
		    
		    
		    self.npp[n+1] = self.sf[n+1] + self.rf[n+1] + self.if_[n+1]
		    self.r2[n+1] = self.M*self.rf[n]
		    self.r1[n+1] = self.rf[n] - self.r2[n]
		    self.p[n+1] = self.if_[n]/self.alpha + self.p[n]
		    self.c[n+1] = self.delta*self.if_[n]
	
		    self.x[n+1] = n
	
		return (self.sf, self.rf, self.if_, self.r1, self.r2, self.npp, self.p, self.c, self.x)

	def plot1(self,show=True, savefig=None):
		"""
			TO DO: Acuerdate de ser consistente con las variables. Tienes que añadirle el self si no
			no van a existir en este metodo.
		"""
		plt.figure(figsize = (15,6))
		plt.plot(self.x,self.rf, label = 'Recovered')
		plt.plot(self.x,self.sf, label = 'Suceptible')
		plt.plot(self.x,self.if_, label = 'Infected')
		plt.plot(self.x,self.npp, label = 'S+I+R', ls = '--')
		plt.legend(loc = 'best', fontsize = 13)
		plt.xlabel('Days', labelpad = 15, fontsize = 15)
		plt.ylabel('Pople number ', labelpad = 15, fontsize = 15)
		plt.title('S-I-R epidemic model', fontsize = 15)

		if isinstance(savefig, str):
			plt.savefig(savefig)

		if show:
			plt.show()

	def plot2(self,show=True, savefig=None):	
		plt.figure(figsize = (20,5))
		plt.subplot(2,2,1)
		plt.plot(self.x,self.c, label= 'Beds ', color = 'crimson')
		plt.legend(loc = 'best', fontsize = 15)
		plt.xlabel('Days', labelpad = 15, fontsize = 15)
		plt.ylabel('Beds number ', labelpad = 15, fontsize = 15)
	
		plt.subplot(2,2,2)
		plt.plot(self.x,self.p, label= 'Positive test (R1)', color = 'k')
		plt.legend(loc = 'best', fontsize = 15)
		plt.xlabel('Days', labelpad = 15, fontsize = 15)
		plt.ylabel('Possitive tests', labelpad = 15, fontsize = 15)
	
	
		plt.subplot(2,2,3)
		plt.plot(self.x,self.r2, label= 'People Death (R2)', color = 'magenta')
		plt.legend(loc = 'best', fontsize = 15)
		plt.xlabel('Days', labelpad = 15, fontsize = 15)
		plt.ylabel('Pople number ', labelpad = 15, fontsize = 15)
	
	
		plt.subplot(2,2,4)
		plt.plot(self.x,self.r1, label= 'People recovered (R1)', color = 'olive')
		plt.legend(loc = 'best', fontsize = 15)
		plt.xlabel('Days', labelpad = 15, fontsize = 15)
	
		if isinstance(savefig, str):
			plt.savefig(savefig)

		if show:
			plt.show()


# MAIN PROGRAM
def main():

	# Crea la instancia EulerRichardson
	data = EulerRichardson(312,0.1,100000,0,1,3)

	# Aplica el metodo para resolver la ED
	data.solve()	
	
	# Grafica los resultados
	data.plot1(show=True, savefig='SIR_Eu_Rich.png')
	data.plot2(show=True, savefig='beds_positive_test_200I.png')


# La siguiente sentencia evita que se ejecute alguna funcion
# cuando el script es importado y no ejecutado desde la terminal
if __name__ == "__main__":
	main()


