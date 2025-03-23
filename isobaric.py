from specific_heat import specific_heats
from entropy import *
from state import State


def isobaric(RS, prop, typ, R, const_cap):

	P = RS.P
	cp, cv, k = specific_heats(RS.T, R, const_cap)

	match typ:
		case 'v':
			v = prop
			T = P * v / R
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_var_cap(RS, T, P, R)
		case 'T':
			T = prop
			v = R * T / P
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_var_cap(RS, T, P, R)
		case _:
			assert(False)

	return State(P, v, T, s)