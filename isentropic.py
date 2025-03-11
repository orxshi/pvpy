from state import State
from specific_heat import specific_heats


def isentropic(RS, prop, typ, R, const_cap):
	'''
	Determines other thermodynamic properties from the given one using isentropic polytropic process.

	For a isentropic polytropic process, PV ** k = constant = C, where, k is the ratio of specific heat capacities.
	The constant is obtained from reference state: C = RS.P * RS.v ** k

	Arguments:
		RS (State): Reference state from which constant of polytropic process is obtained
		prop (float): Numerical value of known property
		typ (str): Type of known property: One of 'P', 'v', 'T'
		R (float): Specific gas constant (K/kg.K)
		const_cap: True if speficic heat capacity to be taken as constant

	Returns:
		State (State): The thermodynamic state.
	'''

	s = RS.s
	cp, cv, k = specific_heats(RS.T, R, const_cap)
	C = RS.P * RS.v ** k

	match typ:
		case 'P':
			P = prop
			v = (C / P) ** (1 / k)
			T = P * v / R
		case 'v':
			v = prop
			P = C / v ** k # wrong! works only for const cap!
			T = P * v / R
		case 'T':
			T = prop
			# cp, cv, k = specific_heats(T, R, const_cap)
			# v = (C / R * T) ** (1 / (k - 1))
			v = RS.v / (T / RS.T) ** (1 / (k - 1))
			P = R * T / v 
		case _:
			assert(False)

	return State(P, v, T, s)