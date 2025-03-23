from state import State
from specific_heat import specific_heats
from entropy import *
from math import exp

def isothermal(RS, prop, typ, R, const_cap):
	'''
	Determines other thermodynamic properties from the given one using isothermal polytropic process.

	For a isothermal polytropic process, PV = constant = C.
	The constant is obtained from reference state: C = RS.P * RS.v

	Arguments:
		RS (State): Reference state from which constant of polytropic process is obtained
		prop (float): Numerical value of known property
		typ (str): Type of known property: One of 'P', 'v', 'T'
		R (float): Specific gas constant (K/kg.K)
		const_cap: True if speficic heat capacity to be taken as constant

	Returns:
		State (State): The thermodynamic state.
	'''

	T = RS.T
	cp, cv, k = specific_heats(RS.T, R, const_cap)
	C = RS.P * RS.v

	match typ:
		case 'P':
			P = prop
			v = C / P
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_var_cap(RS, T, P, R)
		case 'v':
			v = prop
			P = C / v
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_var_cap(RS, T, P, R)
		case 's':
			s = prop
			v = RS.v * exp((s - RS.s) / R)
			P = C / v
		case _:
			assert(False)

	return State(P, v, T, s)