class State:
	'''
	Thermodynamic state. Requires at least two intensive properties. Since entropy is meaningful in differences, entropy can be given as zero.

	Arguments:
		P (float): Absolute pressure (Pa)
		v (float): Specific volume (kg/m3)
		s (float): Entropy (J/kg.K)
		R (float): Specific gas constant (J/kg.K)

	Constructors:
		PT(cls, P, T, R): Makes a State from P (pressure) and T (temperature). Specific gas constant (R) is also required for ideal gas equation.
		
		Pv(cls, P, v, R): Makes a State from P (pressure) and v (specific volume). Specific gas constant (R) is also required for ideal gas equation.

	'''

	def __init__(self, P=None, v=None, T=None, s=None):

		self.P = P
		self.v = v
		self.T = T
		self.s = s

	@classmethod
	def PT(cls, P, T, s, R):

		v = R * T / P
		return State(P, v, T, s)
	
	@classmethod
	def Pv(cls, P, v, s, R):

		T = P * v / R
		return State(P, v, T, s)
	
	def __str__(self):

		return "P: {:.1f}\nv: {:.1f}\nT: {:.1f}\ns: {:.1f}".format(self.P, self.v, self.T, self.s)