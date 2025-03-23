from state import State
from specific_heat import specific_heats
from entropy import *


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
		const_cap: True if specific heat capacity to be taken as constant

	Returns:
		State (State): The thermodynamic state.
	'''

	s = RS.s
	cp, cv, k = specific_heats(RS.T, R, const_cap)
	C = RS.P * RS.v ** k

	match typ:
		case 'P':
			P = prop
			if const_cap:
				T = T_isentropic_Pgiven_const_cap(P, RS, k)
			else:
				T = T_isentropic_Pgiven_var_cap(P, RS, R)
			v = R * T / P
		case 'v':
			v = prop
			if const_cap:
				T = T_isentropic_vgiven_const_cap(v, RS, k)
			else:
				T = T_isentropic_vgiven_var_cap(v, RS, R)
			P = R * T / v
		case 'T':
			T = prop
			if const_cap:
				v = v_isentropic_Tgiven_const_cap(T, RS, k)
			else:
				v = v_isentropic_Tgiven_var_cap(T, RS, R)
			P = R * T / v
		case _:
			assert(False)

	return State(P, v, T, s)