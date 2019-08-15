#INTERPOLATION
def mne_newton(vec_x, vec_y):
    matrix = [(len(vec_x) + 1)*[0] for i in range(len(vec_x))]

    for i in range(len(vec_x)):
        matrix[i][0] = vec_x[i]
        matrix[i][1] = vec_y[i]

    for j in range(2, len(vec_x) + 1):
        for i in range(j-1, len(vec_x)):
            matrix[i][j] = ( matrix[i][j-1] - matrix[i-1][j-1] ) / (matrix[i][0] - matrix[i-j+1][0])

    operators_dif = []
    for i in range(len(vec_x)):
        operators_dif.append(matrix[i][i+1])

    return get_coef_polyn(operators_dif, vec_x)

def get_coef_polyn(consts, vec_x):

	polyn = [consts[0]]
	actual = [1]

	for i in range(1, len(consts)):

		actual.insert(0, 0)
		aux1 = actual[-1]
		for k in reversed(range(len(actual) - 1)):
			aux2 = actual[k]
			actual[k] -= vec_x[i-1] * aux1
			aux1 = aux2

			polyn[k] += consts[i] * actual[k]

		polyn.append(consts[i])

	return polyn

def func(x, y):
	return y

def euler_method(f_xy, x0, y0, h, num_points):

	x = [x0]
	y = [y0]

	for i in range(num_points - 1):
		y.append(y[-1] + h * f_xy(x[-1], y[-1]) )
		x.append(x[-1] + h)

	#Return a vector with coefficients of a interpolator polynomial
	return mne_newton(x, y)



import matplotlib.pyplot as plt
import numpy as np
import math

x = np.arange(-0.1, 0.6, 0.01)
y1 = np.poly1d(euler_method(func, 0, 1, 0.1, 5)[::-1])
plt.plot(x, y1(x))
plt.plot(x, np.exp(x))

plt.show()
