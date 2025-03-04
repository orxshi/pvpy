from math import log as ln
from std_entropy import snot


def entropy_const_cap(RS, T, v, cv, R):
	return RS.s + cv * ln(T/RS.T) + R * ln(v/RS.v)

def entropy_var_cap(RS, T, v, R):
	return snot(T) - snot(RS.T) + R * ln(v/RS.v)