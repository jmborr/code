from mantid.api import *
import numpy as np

# Class definition
class QHO_Q_function(IFunction1D):
 	# Q dependence of the scattering from a Quantum Harmonic Oscillator for Nth level

	def init(self):
		self.declareParameter("Amplitude",0.01)
		self.declareParameter("Constant",0.0)
		self.declareParameter("W",0.02)
		self.declareAttribute("Order",0)
		
	def setAttributeValue(self, name, value):
		if  name== "Order":
			self._N = value
			
	def function1D(self,xvals):
		# Load parameters
		A=self.getParameterValue("Amplitude")	
		C = self.getParameterValue("Constant")
		w=self.getParameterValue("W")	
		N = self._N

		return A*xvals**(2*N)*np.exp(-w*xvals**2) + C 

# Register algorithm with Mantid
FunctionFactory.subscribe(QHO_Q_function)

