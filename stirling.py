from state import State
from isochoric import isochoric
from isothermal import isothermal
from specific_heat import specific_heats, get_temperature_by_integrating_cv_dT, get_heat_by_integrating_cv_dT
from myarray import array
from print import *
from work import work_qr_qs
from pvts import *
from conversion import celsius
from copy import copy


def stirling_Tmax_qs(Tmax, qs, GS, const_cap):

	states = []

	R = 287 # J/kg.K
	cp, cv, k = specific_heats(GS.T, R, const_cap)

	RS = GS

	# isochoric expansion
	for T in array(GS.T, Tmax):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	s3 = RS.s - qs / Tmax

	# isothermal heat addition
	for s in array(RS.s, s3):
		S = isothermal(RS, s, 's', R, const_cap)
		states.append(S)

	RS = states[-1]

	# compression ratio
	cr = GS.v / RS.v
	
	# isochoric compression
	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	# get heat rejected qr
	qr = abs(RS.T * (RS.s - GS.s))

	# isothermal heat rejection
	for v in array(RS.v, GS.v):
		S = isothermal(RS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]
	
	w = work_qr_qs(qr, qs)

	printres('Stirling', w, qs, qr, cr)

	return states


def stirling(cr, qs, GS, const_cap):
	'''
	Stirling cycle.
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

	# isothermal heat rejection
	for v in array(GS.v, vL):
		S = isothermal(GS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]

	# get heat rejected qr
	qr = abs(RS.T * (RS.s - GS.s))

	# get Tmax
	Tmax = qs / (GS.s - RS.s)

	# isochoric compression
	for T in array(RS.T, Tmax):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)
	
	RS = states[-1]

	# isothermal heat addition
	for v in array(vL, GS.v):
		S = isothermal(RS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]

	# isochoric expansion
	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)
	
	w = work_qr_qs(qr, qs)

	printres('Stirling', w, qs, qr, cr)

	return states


# GS = State.PT(101325, 300, 0, 287)
# states = stirling(9.890421999965, 96000, GS, True)
# states = stirling_Tmax_qs(GS.T + 300, 96000, GS, True)


# pv([states], ['Stirling'])
# ts([states], ['Stirling'])
