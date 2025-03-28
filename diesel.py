from state import State
from isentropic import isentropic
from isochoric import isochoric
from specific_heat import specific_heats
from isobaric import isobaric
from myarray import array
from work import work
from print import *
from pvts import *
from copy import copy

def diesel_Pmax_Tmax_qr(Pmax, Tmax, qr, const_cap):

	states = []

	R = 287 # J/kg.K

	S3 = State.PT(Pmax, Tmax, 0, R)
	RS = S3

	cp, cv, k = specific_heats(RS.T, R, const_cap)

	# isentropic expansion
	for v in array(RS.v, RS.v + 1.0):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]
	S4 = copy(RS)

	T1 = RS.T - qr / cv

	# constant-volume heat rejection
	for T in array(RS.T, T1):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	# constant-volume heat rejection
	for T in array(RS.T, T1):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

	# isentropic compression
	for P in array(RS.P, S3.P):
		S = isentropic(RS, P, 'P', R, const_cap)
		states.append(S)

	RS = states[-1]

	qs = cp * (S3.T - RS.T)

	cr = S4.v / RS.v

	# isobaric heat addition
	for T in array(RS.T, S3.T):
		S = isobaric(RS, T, 'T', R, const_cap)
		states.append(S)

	w = work(states)

	printres('Diesel', w, qs, qr, cr)

	return states
	


def diesel(cr, qs, GS, const_cap):
	'''
	Diesel cycle.
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

	# isentropic compression
	for v in array(RS.v, vL):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)
	
	RS = states[-1]
	
    # intermediate temperature at the end of constant-pressure heat addition
	if const_cap is True:
		T3 = RS.T + qs / cp
	else:
		pass
	
    # isobaric heat addition
	for T in array(RS.T, T3):
		S = isobaric(RS, T, 'T', R, const_cap)
		states.append(S)

	RS = states[-1]

    # isentropic expansion
	for v in array(RS.v, GS.v):
		S = isentropic(RS, v, 'v', R, const_cap)
		states.append(S)

	RS = states[-1]

    # isochoric heat rejection
	for T in array(RS.T, GS.T):
		S = isochoric(RS, T, 'T', R, const_cap)
		states.append(S)
		
    # get heat rejected qr
	if const_cap is True:
		qr = cv * (RS.T - GS.T)
	else:
		pass

	w = work(states)

	printres('Diesel', w, qs, qr, cr)
	
	return states


# GS = State.PT(101325, 300, 0, 287)
# states = diesel(6, 1170000, GS, True)

# states = diesel_Pmax_Tmax_qr(200000, 600, 60000, True)

# pv([states], ['Diesel'])
# ts([states], ['Diesel'])