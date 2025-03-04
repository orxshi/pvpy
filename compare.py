from otto import otto
from diesel import diesel
from carnot import carnot
from pvts import *
from state import State


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


# same_cr_same_heat_input_const_cap()
same_cr_same_heat_input_const_cap_carnot()