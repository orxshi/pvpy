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


def isentropic(RS, prop, typ, k, R):

	s = RS.s
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


def pv(states, label):
	plt.plot([s.v for s in states], [s.P for s in states], label = label)


def work(states):
	w = 0

	for i in range(1, len(states)):
		w += states[i].P * (states[i].v - states[i-1].v)

	return w


def otto(cr, QS):

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
	T3 = RS.T + QS / cv

	for T in array(RS.T, T3):
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
	
	print('Work from Otto = {:.0f} kJ'.format(work(states)/1000))

	return states


def diesel(cr, QS):

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
	T3 = RS.T + QS / cp

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
	
	print('Work from Diesel = {:.0f} kJ'.format(work(states)/1000))

	return states


def main():
	states_otto = otto(6, 500000)
	states_diesel = diesel(6, 500000)
	pv(states_otto, 'Otto')
	pv(states_diesel, 'Diesel')
	plt.show()
	

if __name__ == '__main__':
	main()












