from state import State
from isentropic import isentropic
from isochoric import isochoric
from specific_heat import specific_heats
from isobaric import isobaric
from myarray import array
from work import work
from print import *
from pvts import *


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

	printres('Diesel', w, qs, qr)
	
	return states


# GS = State.PT(101325, 300, 0, 287)
# states = diesel(6, 1170000, GS, True)

# pv([states], ['Diesel'])
# ts([states], ['Diesel'])