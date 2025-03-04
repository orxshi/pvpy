def work(states):
	w = 0

	for i in range(1, len(states)):
		w += states[i].P * (states[i].v - states[i-1].v)

	return w