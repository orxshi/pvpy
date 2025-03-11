from specific_heat import specific_heats
from entropy import *
from state import State


def isochoric(RS, prop, typ, R, const_cap):

	v = RS.v
	cp, cv, k = specific_heats(RS.T, R, const_cap)

	match typ:
		case 'P':
			P = prop
			T = P * v / R
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case 'T':
			T = prop
			P = R * T / v
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case _:
			assert(False)

	return State(P, v, T, s)