import numpy as np
from copy import copy

def specific_heats(T, R, const_cap):
	'''
	Returns specific heat capacity at
	
	Arguments:
        T (float): Temperature in Kelvin
		R (float): Specific gas constant (J/kg.K)
		const_cap (bool): True if specific heat capacity is constant.
		
	Returns:
        cp (float), cv (float), k (float):
		Specific heat capacity at constant pressure (J/kg.K),
		specific heat capacity at constant volume (J/kg.K),
		specific heat capacity ratio (cp/cv)
	'''

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


def get_temperature_by_integrating_cv_dT(RS, R, q):
	'''
	Finds temperature from given heat input by integrating du = cv * dT.
	When there is no work done by/to control mass, dq = du.

	Returns:
		T (float): The temperature which corresponds to given amount of heat if the temperature is less than Tmax=5000 K. Otherwise, returns None.

	Nomenclature:
		qcur: Current amount of heat
	'''
	qcur = 0
	prevT = RS.T
	Tmax = 5000
	for T in np.arange(RS.T, Tmax, 1000):
		cp, cv, k = specific_heats(RS.T, R, False)
		qcur += cv * (T - prevT)
		prevT = copy(T)
		if qcur >= q:
			return T
	
	return None

def get_heat_by_integrating_cv_dT(Ta, Tb, R):
	q = 0
	prevT = Ta
	for T in np.arange(Ta, Tb, 1000):
		cp, cv, k = specific_heats(T, R, False)
		q += cv * (T - prevT)
		prevT = copy(T)

	return q

