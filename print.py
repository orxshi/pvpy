def printres(cycle, w, qs, qr):

	err_w = abs(w - (qs - qr)) / w * 100
	# ef = w / qs * 100
	ef = (1 - qr / qs) * 100

	print(cycle)
	print('----------------')
	print('w = {:.0f} kJ/kg'.format(w / 1000))
	print('qs = {:.0f} kJ/kg'.format(qs / 1000))
	print('qr = {:.0f} kJ/kg'.format(qr / 1000))
	print('ef = {:.0f}%'.format(ef))
	print('err_w = {:.2f}%'.format(err_w))
	print('----------------')