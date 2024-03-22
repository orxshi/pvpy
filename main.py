from numpy import linspace as linspace
from math import log as ln
import matplotlib.pyplot as plt

class State:

	def __init__(self, P, v, s, R):

		self.P = P
		self.v = v
		self.T = P * v / R
		self.s = s
	
	def __str__(self):

		return "P: {:.1f}\nv: {:.1f}\nT: {:.1f}\ns: {:.1f}".format(self.P, self.v, self.T, self.s)


def entropy(RS, T, v, cv, R):
	return RS.s + cv * ln(T/RS.T) + R * ln(v/RS.v)


def specific_heats(T, R):

	# T must be in Kelvin.

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	H = T / 1000

	cp = C0 + C1 * H + C2 * H ** 2 + C3 * H ** 3
	cv = cp - R
	k = cp / cv

	return cp, cv, k


def isentropic(RS, prop, typ, k, R):

	s = RS.s
	cp, cv, k = specific_heat_p(RS.T)
	C = RS.P * RS.v ** k

	match typ:
		case 'P':
			P = prop
			v = (C / P) ** (1 / k)
			T = P * v / R
		case 'v':
			v = prop
			P = C / v ** k
			T = P * v / R
		case 'T':
			T = prop
			cp, cv, k = specific_heat_p(T)
			v = (C / R * T) ** (1 / (k - 1))
			P = R * T / v 
		case _:
			assert(False)

	return State(P, v, s, R)


def isobaric(RS, prop, typ, cv, R):

	P = RS.P

	match typ:
		case 'v':
			v = prop
			T = P * v / R
			s = entropy(RS, T, v, cv, R)
		case 'T':
			T = prop
			v = R * T / P
			s = entropy(RS, T, v, cv, R)
		case _:
			assert(False)

	return State(P, v, s, R)


def isochoric(RS, prop, typ, cv, R):

	v = RS.v

	match typ:
		case 'P':
			P = prop
			T = P * v / R
			s = entropy(RS, T, v, cv, R)
		case 'T':
			T = prop
			P = R * T / v
			s = entropy(RS, T, v, cv, R)
		case _:
			assert(False)

	return State(P, v, s, R)


def array(a, b):
	return linspace(a, b, 4320)


def pv(states, labels):

	plt.xlabel(r"$\nu$ (m$^3/$kg)")
	plt.ylabel("P (kPa)")
	plt.title(r"P-$\nu$ diagram")

	for i in range(len(states)):
		plt.plot([s.v for s in states[i]], [s.P/1000 for s in states[i]], label = labels[i])

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


def otto(cr, qs):

	states = []

	R = 287

	GS = State(101325, 0.6, 0, R)
	RS = GS
	vL = RS.v / cr

	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', k, R)
		states.append(S)
	
	RS = states[-1]
	cp, cv, k = specific_heat_p(RS.T)
	T3 = RS.T + qs / cv

	for T in array(RS.T, T3):
		cp, cv, k = specific_heat_p(T)
		S = isochoric(RS, T, 'T', cv, R)
		states.append(S)

	RS = states[-1]

	for v in array(vL, GS.v):
		S = isentropic(RS, v, 'v', k, R)
		states.append(S)

	RS = states[-1]

	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', cv, R)
		states.append(S)
	
	qr = cv * (RS.T - GS.T)

	w = work(states)

	printres('Otto', w, qs, qr)

	return states


def diesel(cr, qs):

	states = []

	cp = 1004
	R = 287
	cv = cp - R
	k = cp / cv

	GS = State(101325, 0.6, 0, R)
	RS = GS
	vL = RS.v / cr

	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', k, R)
		states.append(S)
	
	RS = states[-1]
	T3 = RS.T + qs / cp

	for T in array(RS.T, T3):
		S = isobaric(RS, T, 'T', cv, R)
		states.append(S)

	RS = states[-1]

	for v in array(RS.v, GS.v):
		S = isentropic(RS, v, 'v', k, R)
		states.append(S)

	RS = states[-1]

	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', cv, R)
		states.append(S)

	qr = cv * (RS.T - GS.T)

	w = work(states)

	printres('Diesel', w, qs, qr)
	
	return states


def same_cr_same_heat_input():

	cr = 6
	qs = 500000

	print('Results for given cr = {:.0f} and qs = {:.0f} kJ/kg'.format(cr, qs/1000))
	print('----------------')

	states_otto = otto(cr, qs)
	states_diesel = diesel(cr, qs)

	pv([states_otto, states_diesel], ['Otto', 'Diesel'])


def main():

	same_cr_same_heat_input()
	

if __name__ == '__main__':
	main()













