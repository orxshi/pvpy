from specific_heat import specific_heats
from entropy import *
from state import State


def isobaric(RS, prop, typ, R, const_cap):

	P = RS.P

	match typ:
		case 'v':
			v = prop
			T = P * v / R
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case 'T':
			T = prop
			v = R * T / P
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case _:
			assert(False)

	return State(P, v, T, s)