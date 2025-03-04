from state import State
from isentropic import isentropic
from isothermal import isothermal
from specific_heat import specific_heats, get_temperature_by_integrating_cv_dT, get_heat_by_integrating_cv_dT
from myarray import array
from print import *
from work import work
from pvts import *
from conversion import celsius
from copy import copy


def carnot(cr, qs, GS, const_cap):
	'''
	Carnot cycle.
	TODO: const_cap = False does not work since isentropic.py not ready for variable specific heat capacity. 

	Arguments:
		cr (float): Compression ratio
		qs (float): Specific heat supplied / input (J/kg)
		GS: Ground state which corresponds to beginning of compression stroke
		const_cap: True if specific heat capacity is to be constant

	Returns:
		states ([State]): List of states in the cycle.

	Nomenclature:
		RS: Reference state which is the last state of each process
		vL: Low specific volume which corresponds to end of compression stroke
		qr: Rejected/released heat
	'''

	states = []

	R = 287 # J/kg.K
	cp, cv, k = specific_heats(GS.T, R, const_cap)

	RS = GS
	vL = RS.v / cr


	# get state 2 with iteration
	error = 10000
	tol = 1000
	vc = vL
	while abs(error) > tol:
		# isothermal heat rejection
		vc += vL / 100
		for v in array(GS.v, vc):
			S = isothermal(GS, v, 'v', R, const_cap)
		
		RS = S
		MS = copy(S)

		# isentropic compression
		for v in array(RS.v, vL):
			S = isentropic(RS, v, 'v', R, const_cap)

		RS = S

		qs_calc = RS.T * (GS.s - RS.s)
		error = qs_calc - qs
		# print(error)

	RS = MS

	# print(MS.v, MS.T, MS.P, MS.s)

	GS.s -= RS.s
	RS.s -= RS.s








	











		




	# isothermal heat rejection
	for s in array(GS.s, RS.s):
		S = isothermal(GS, s, 's', R, const_cap)
		states.append(S)

	RS = states[-1]

	# get heat rejected qr
	qr = abs(RS.T * (RS.s - GS.s))

	# isentropic compression
	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)
	
	RS = states[-1]

	# isothermal heat addition
	for s in array(RS.s, GS.s):
		S = isothermal(RS, s, 's', R, const_cap)
		states.append(S)

	RS = states[-1]

	# isentropic expansion
	for v in array(RS.v, GS.v):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)
	
	w = work(states)

	printres('Carnot', w, qs, qr)

	return states


# Ganesan worked-out 2.7
# GS = State.PT(101325, 300, 0, 287)
# states = carnot(6, 117000, GS, True)
# carnot(6, 117000, GS, True)

# pv([states], ['Carnot'])
# ts([states], ['Carnot'])