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
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_ana_var_cap(T, P, RS, R)
				# s = entropy_var_cap(RS, T, P, R)
		case 'T':
			T = prop
			P = R * T / v
			s = entropy_const_cap(RS, T, v, cv, R)
			if const_cap == False:
				s = entropy_ana_var_cap(T, P, RS, R)
				# s = entropy_var_cap(RS, T, P, R)
		case _:
			assert(False)

	return State(P, v, T, s)