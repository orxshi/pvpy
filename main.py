from math import log as ln
import matplotlib.pyplot as plt
from std_entropy import snot
import numpy as np
import copy





class State:

	def __init__(self, P, v, s, R):

		self.P = P
		self.v = v
		self.T = P * v / R
		self.s = s
	
	def __str__(self):

		return "P: {:.1f}\nv: {:.1f}\nT: {:.1f}\ns: {:.1f}".format(self.P, self.v, self.T, self.s)


def entropy_const_cap(RS, T, v, cv, R):
	return RS.s + cv * ln(T/RS.T) + R * ln(v/RS.v)

def entropy_var_cap(RS, T, v, R):
	return snot(T) - snot(RS.T) + R * ln(v/RS.v)


def specific_heats(T, R, const_cap):

	# T must be in Kelvin.

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	H = T / 1000

	cp = C0 + C1 * H + C2 * H ** 2 + C3 * H ** 3
	cp *= 1000

	if const_cap is True:
		cp = 1005

	cv = cp - R
	k = cp / cv

	return cp, cv, k


def isentropic(RS, prop, typ, R, const_cap):

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
			cp, cv, k = specific_heats(T, R, const_cap)
			v = (C / R * T) ** (1 / (k - 1))
			P = R * T / v 
		case _:
			assert(False)

	return State(P, v, s, R)


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

	return State(P, v, s, R)


def isochoric(RS, prop, typ, R, const_cap):

	v = RS.v

	match typ:
		case 'P':
			P = prop
			T = P * v / R
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case 'T':
			T = prop
			P = R * T / v
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			if const_cap is True:
				s = entropy_const_cap(RS, T, v, cv, R)
			else:
				s = entropy_var_cap(RS, T, v, R)
		case _:
			assert(False)

	return State(P, v, s, R)


def array(a, b):
	return np.linspace(a, b, 4320)


def pv(states, labels):

	plt.xlabel(r"$\nu$ (m$^3/$kg)")
	plt.ylabel("P (kPa)")
	plt.title(r"P-$\nu$ diagram")

	for i in range(len(states)):
		plt.plot([s.v for s in states[i]], [s.P/1000 for s in states[i]], label = labels[i])

	plt.legend(loc='best')
	plt.show()


def ts(states, labels):

	plt.xlabel(r"s (J/kgK)")
	plt.ylabel("T (K)")
	plt.title(r"T-s diagram")

	for i in range(len(states)):
		plt.plot([s.s for s in states[i]], [s.T for s in states[i]], label = labels[i])

	plt.legend(loc='best')
	plt.show()


def work(states):
	w = 0

	for i in range(1, len(states)):
		w += states[i].P * (states[i].v - states[i-1].v)

	return w


def printres(cycle, w, qs, qr):

	err_w = abs(w - (qs - qr)) / w * 100
	ef = w / qs * 100

	print(cycle)
	print('----------------')
	print('w = {:.0f} kJ/kg'.format(w / 1000))
	print('qs = {:.0f} kJ/kg'.format(qs / 1000))
	print('qr = {:.0f} kJ/kg'.format(qr / 1000))
	print('ef = {:.0f}%'.format(ef))
	print('err_w = {:.2f}%'.format(err_w))
	print('----------------')


def otto(cr, qs, const_cap):

	states = []

	R = 287

	GS = State(101325, 0.6, 0, R)
	RS = GS
	vL = RS.v / cr

	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)
	
	RS = states[-1]

	if const_cap is True:
		cp, cv, k = specific_heats(RS.T, R, const_cap)
		T3 = RS.T + qs / cv
	else:
		qsn = 0
		prevT = RS.T
		for T in np.arange(RS.T, 999999999, 10):
			S = isochoric(RS, T, 'T', R, const_cap)
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			qsn += cv * (T - prevT)
			prevT = copy.copy(T)
			if qsn >= qs:
				T3 = T
				break
	

	print('T3: ', T3)
	
	for T in array(RS.T, T3):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	for v in array(vL, GS.v):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]

	qr = 0
	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)
		if const_cap is False:
			cp, cv, k = specific_heats(RS.T, R, const_cap)
			qr -= cv * (S.T - states[-2].T)
	
	if const_cap is True:
		qr = cv * (RS.T - GS.T)

	w = work(states)

	printres('Otto', w, qs, qr)

	return states


def diesel(cr, qs, const_cap):

	states = []

	R = 287

	GS = State(101325, 0.6, 0, R)
	RS = GS
	vL = RS.v / cr

	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', R)
		states.append(S)
	
	RS = states[-1]
	T3 = RS.T + qs / cp

	for T in array(RS.T, T3):
		S = isobaric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	for v in array(RS.v, GS.v):
		S = isentropic(RS, v, 'v', R)
		states.append(S)

	RS = states[-1]

	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', cv, R)
		states.append(S)

	qr = cv * (RS.T - GS.T)

	w = work(states)

	printres('Diesel', w, qs, qr)
	
	return states


def same_cr_same_heat_input_const_cap():

	cr = 6
	qs = 500000

	print('Results for given cr = {:.0f} and qs = {:.0f} kJ/kg'.format(cr, qs/1000))
	print('----------------')

	states_otto = otto(cr, qs, const_cap = True)
	states_diesel = diesel(cr, qs, const_cap = True)

	pv([states_otto, states_diesel], ['Otto', 'Diesel'])
	ts([states_otto, states_diesel], ['Otto', 'Diesel'])


def same_cr_same_heat_input():

	cr = 6
	qs = 500000

	print('Results for given cr = {:.0f} and qs = {:.0f} kJ/kg'.format(cr, qs/1000))
	print('----------------')

	states_const = otto(cr, qs, const_cap = True)
	states_var = otto(cr, qs, const_cap = False)
	#states_diesel = diesel(cr, qs, const_cap = True)

	pv([states_const, states_var], ['const', 'var'])
	#ts([states_otto, states_diesel], ['Otto', 'Diesel'])


def main():

	same_cr_same_heat_input()
	

if __name__ == '__main__':
	main()
