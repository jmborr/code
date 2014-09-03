from mantid.api import *
import numpy as np
import scipy.special as sps


def ICGauss(a,b,R,s,t) :
    # Based on Jack Carpenter notes from "Elements of Slow Neutron Scattering"
    # Modified to take Gaussian sigma directly using recursive formula (DLA)

    prefactor = a/2.0/np.sqrt(2.0*np.pi)/s
    exp_a = np.exp(-a*t)
    c1 = (1.0 - R)*a**2
    d = a - b    
    c2 = 2.0*R*(a**2*b)/d**3

    prefactor_gn = np.sqrt(np.pi/2.0)*s
    ga0 = prefactor_gn*np.exp(0.5*a**2*s**2)*sps.erfc((a*s - t/s)/np.sqrt(2.0))
    gb0 = prefactor_gn*np.exp(0.5*b**2*s**2)*sps.erfc((b*s - t/s)/np.sqrt(2.0))
    ga1 = s**2*(-a*ga0 + np.exp(a*t - t**2/2.0/s**2))
    ga2 = s**2*(ga0 - a*ga1 - t*np.exp(a*t - t**2/2.0/s**2))    

    term1 = (c1*t**2 - c2*(1.0 + d*t + 0.5*d**2*t**2))*exp_a*ga0
    term2 = c2*np.exp(-b*t)*gb0
    term3 = (2*c1*t - c2*d*(1.0+d*t))*exp_a*ga1
    term4 = (c1 - c2*d**2/2.0)*exp_a*ga2

    return prefactor*(term1 + term2 + term3 + term4)


# Class definition
class ICxGmonitor(IFunction1D):
 
	def init(self):
		self.declareParameter("Amplitude",1.0)
		self.declareParameter("Tpeak",0.0)
		self.declareParameter("Alpha",0.2)
		self.declareParameter("Beta",0.04)
		self.declareParameter("R",0.6)
		self.declareParameter("Sigma",10.0)
		
	def function1D(self,xvals):
		# Load parameters
		a0=self.getParameterValue("Amplitude")	
		tp = self.getParameterValue("Tpeak")
		a=self.getParameterValue("Alpha")	
		b = self.getParameterValue("Beta")
		R=self.getParameterValue("R")	
		s = self.getParameterValue("Sigma")

		tmean = 3.0/a + R/b

		return a0*ICGauss(a,b,R,s,tp+tmean-xvals) 

# Register algorithm with Mantid
FunctionFactory.subscribe(ICxGmonitor)

