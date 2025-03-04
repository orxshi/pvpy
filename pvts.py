import matplotlib.pyplot as plt


def pv(states, labels):

	plt.xlabel(r"$\nu$ (m$^3/$kg)")
	plt.ylabel("P (kPa)")
	plt.title(r"P-$\nu$ diagram")

	for i in range(len(states)):
		plt.plot([s.v for s in states[i]], [s.P/1000 for s in states[i]], label = labels[i])

	plt.legend(loc='best')
	plt.show()


def ts(states, labels):

	plt.xlabel(r"s (J/kgK)")
	plt.ylabel("T (K)")
	plt.title(r"T-s diagram")

	for i in range(len(states)):
		plt.plot([s.s for s in states[i]], [s.T for s in states[i]], label = labels[i])

	plt.legend(loc='best')
	plt.show()