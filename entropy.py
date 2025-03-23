from math import log as ln
from math import exp
from std_entropy import snot
from copy import copy


def entropy_const_cap(RS, T, v, cv, R):
	return RS.s + cv * ln(T/RS.T) + R * ln(v/RS.v)

def T_isentropic_Pgiven_const_cap(P, RS, k):
	return RS.T * (P / RS.P) ** (1 - 1 / k)

def T_isentropic_vgiven_const_cap(v, RS, k):
	return RS.T * (RS.v / v) ** (k - 1)

def v_isentropic_Tgiven_const_cap(T, RS, k):
	return RS.v * (RS.T / T) ** (1 / (k - 1))

# def entropy_var_cap(RS, T, v, R):
	# return snot(T) - snot(RS.T) + R * ln(v/RS.v)

def entropy_var_cap(RS, T, P, R):
	return snot(T) - snot(RS.T) - R * ln(P/RS.P)

def entropy_ana_var_cap(T, P, RS, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000
	
	delta_s = C0 * ln(T / RS.T)\
		+ C1 * (T ** 1 - RS.T ** 1) / 1e3\
		+ C2 * (T ** 2 - RS.T ** 2) / 2e6\
		+ C3 * (T ** 3 - RS.T ** 3) / 3e9\
		- R * ln(P / RS.P)
	
	# delta_s = 1.005\
	# 		- R * ln(P / RS.P)
	
	return RS.s + delta_s

def specific_volume_isentropic_var_cap(T, RS, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000

	A = (C0 - R) * ln(T / RS.T) + C1 * (T - RS.T) / 1000 + C2 * (T ** 2 - RS.T ** 2) / 2e6 + C3 * (T ** 3 - RS.T ** 3) / 3e9

	v = exp(-A / R) * RS.v

	return v

def T_isentropic_vgiven_var_cap(v, RS, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000

	T = RS.T # initial guess
	error = float('inf')
	tolerance = 1

	while error > tolerance:

		f = (C0 - R) * ln(T / RS.T)\
		+ C1 * (T ** 1 - RS.T ** 1) / 1e3\
		+ C2 * (T ** 2 - RS.T ** 2) / 2e6\
		+ C3 * (T ** 3 - RS.T ** 3) / 3e9\
		+ R * ln(v / RS.v)
		
		fp = (C0 - R) / T\
		+ C1 * T ** 0 / 1e3\
		+ C2 * T ** 1 / 1e6\
		+ C3 * T ** 2 / 1e9
		
		T_old = copy(T)
		T = T - f / fp

		error = abs(T - T_old)
		# print(error, T_old, T, RS.T, f, fp, v, RS.v)

	return T

def T_isentropic_Pgiven_var_cap(P, RS, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000

	T = RS.T # initial guess
	error = float('inf')
	tolerance = 1

	while error > tolerance:

		f = C0 * ln(T / RS.T)\
		+ C1 * (T ** 1 - RS.T ** 1) / 1e3\
		+ C2 * (T ** 2 - RS.T ** 2) / 2e6\
		+ C3 * (T ** 3 - RS.T ** 3) / 3e9\
		- R * ln(P / RS.P)
		
		fp = C0 / T\
		+ C1 * T ** 0 / 1e3\
		+ C2 * T ** 1 / 1e6\
		+ C3 * T ** 2 / 1e9
		
		T_old = T
		T = T - f / fp

		error = abs(T - T_old)

	return T
	

def v_isentropic_Tgiven_var_cap(T, RS, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000

	v = (R - C0) * ln(T / RS.T) / R\
	  - C1 * (T ** 1 - RS.T ** 1) / 1e3 / R\
	  - C2 * (T ** 2 - RS.T ** 2) / 2e6 / R\
	  - C3 * (T ** 3 - RS.T ** 3) / 3e9 / R\
	  - ln(RS.v)
	
	return exp(v)
		
