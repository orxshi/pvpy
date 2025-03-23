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
	Finds temperature from given heat input by integrating du = cv * dT. cv is the specific heat capacity which depends on temperature.
	When there is no work done by/to control mass, dq = du. That is, q = int cv(T) * dT.

	Arguments:
		RS (float): Reference state
		R (float): Specific gas constant
		q (float): Heat transfer during the process

	Returns:
		T (float): The temperature which corresponds to given amount of heat.
	'''

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000
	
	T = RS.T # initial guess
	error = float('inf')
	tolerance = 1

	# q must be physical otherwise the function will not have root so NR will not converge. Use the snippet below to see te root.
	# Ts = np.linspace(RS.T, 3000, 100)
	# f = []
	# import matplotlib.pyplot as plt
	# for T in Ts:
	# 	ff = (C0 - R) * (T - RS.T)\
	# 		+ C1 * (T ** 2 - RS.T ** 2) / 2e3\
	# 		+ C2 * (T ** 3 - RS.T ** 3) / 3e6\
	# 		+ C3 * (T ** 4 - RS.T ** 4) / 4e9\
	# 		- q
	# 	f.append(ff)
	# plt.plot(Ts, f)
	# plt.axhline(0)
	# plt.show()

	while error > tolerance:

		f = (C0 - R) * (T - RS.T)\
		+ C1 * (T ** 2 - RS.T ** 2) / 2e3\
		+ C2 * (T ** 3 - RS.T ** 3) / 3e6\
		+ C3 * (T ** 4 - RS.T ** 4) / 4e9\
		- q

		fp = C0 - R\
		+ C1 * T / 1e3\
		+ C2 * T ** 2 / 1e6\
		+ C3 * T ** 3 / 1e9

		T_old = copy(T)
		T = T - f / fp

		error = abs(T - T_old)

	return T

def get_heat_by_integrating_cv_dT(Ta, Tb, R):

	C0 = 1.05
	C1 = -0.365
	C2 = 0.85
	C3 = -0.39

	C0 *= 1000
	C1 *= 1000
	C2 *= 1000
	C3 *= 1000

	q = (C0 - R) * (Tb - Ta)\
		+ C1 * (Tb ** 2 - Ta ** 2) / 2e3\
		+ C2 * (Tb ** 3 - Ta ** 3) / 3e6\
		+ C3 * (Tb ** 4 - Ta ** 4) / 4e9

	return q

