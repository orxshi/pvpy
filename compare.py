from otto import otto
from diesel import diesel
from carnot import carnot_Tmax_qs
from stirling import stirling_Tmax_qs
from pvts import *
from state import State


def same_Tmax_qs_const_cap():

	cr = 6
	qs = 500000

	GS = State.PT(101325, 300, 0, 287)

	states_stirling = stirling_Tmax_qs(GS.T + 300, 96000, GS, True)
	states_carnot   = carnot_Tmax_qs  (GS.T + 300, 96000, GS, True)

	pv([states_carnot, states_stirling], ['Carnot', 'Stirling'])
	ts([states_carnot, states_stirling], ['Carnot', 'Stirling'])


def same_cr_same_heat_input_const_cap():

	cr = 6
	qs = 500000

	print('Results for given cr = {:.0f} and qs = {:.0f} kJ/kg'.format(cr, qs/1000))
	print('----------------')

	GS = State.PT(101325, 300, 0, 287)

	states_otto = otto(cr, qs, GS, const_cap = True)
	states_diesel = diesel(cr, qs, GS, const_cap = True)

	pv([states_otto, states_diesel], ['Otto', 'Diesel'])
	ts([states_otto, states_diesel], ['Otto', 'Diesel'])


def same_cr_same_heat_input_const_cap_carnot():

	cr = 6
	qs = 117000

	print('Results for given cr = {:.0f} and qs = {:.0f} kJ/kg'.format(cr, qs/1000))
	print('----------------')

	GS = State.PT(101325, 300, 0, 287)

	states_otto = otto(cr, qs, GS, const_cap = True)
	states_carnot = carnot(cr, qs, GS, const_cap = True)

	pv([states_otto, states_carnot], ['Otto', 'Carnot'])
	ts([states_otto, states_carnot], ['Otto', 'Carnot'])


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


same_cr_same_heat_input_const_cap()
# same_cr_same_heat_input_const_cap_carnot()
# same_Tmax_qs_const_cap()