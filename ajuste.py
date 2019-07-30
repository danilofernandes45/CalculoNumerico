
from math import sqrt, log, exp

def row_criteria(matrix):

	alpha = 0
	for i in range(len(matrix)):

		x = 0

		for j in range(i):
			x += abs(matrix[i][j])

		for j in range(i + 1, len(matrix)):
			x += abs(matrix[i][j])

		x /= abs(matrix[i][i])

		alpha = max(alpha, x)

	if(alpha < 1):
		return True

	return False

def col_criteria(matrix):

	alpha = 0
	for j in range(len(matrix)):

		x = 0

		for i in range(j):
			x += abs(matrix[i][j])

		for i in range(j + 1, len(matrix)):
			x += abs(matrix[i][j])

		x /= abs(matrix[i][i])

		alpha = max(alpha, x)

	if(alpha < 1):
		return True

	return False


def sassenfeld_criteria(matrix):

	alpha = [0] * len(matrix)
	for i in range(len(matrix)):

		for j in range(i):
			alpha[i] += alpha[j] * abs(matrix[i][j])

		for j in range(i + 1, len(matrix)):
			alpha[i] += abs(matrix[i][j])

		alpha[i] /= abs(matrix[i][i])

	if(max(alpha) < 1):
		return True

	return False


def convergence_test(matrix):

	if( row_criteria(matrix) or col_criteria(matrix) or sassenfeld_criteria(matrix) ):
		print("Convergência garantida para os métodos iterativos")
		return True

	print("Convergência não garantida para os métodos iterativos")
	return False

def dist(x1, x2):
	sum = 0
	for i in range(len(x1)):
		sum += (x1[i] - x2[i])**2

	return sqrt(sum)

def nonzero_diagonal(matrix):

	for i in range(len(matrix)):
		if( matrix[i][i] == 0 ):
			for j in range(i, len(matrix)):
				if( matrix[j][i] != 0 ):
					aux = matrix[i].copy()
					matrix[i] = matrix[j].copy()
					matrix[j] = aux.copy()

				if( matrix[i][i] == 0):
					for j in range(i):
						if( matrix[j][i] != 0 and matrix[i][j] != 0):
							aux = matrix[i].copy()
							matrix[i] = matrix[j].copy()
							matrix[j] = aux.copy()

	return matrix

def gauss_seidel(matrix, vector, eps, kmax):

	matrix = nonzero_diagonal(matrix)

	x_1 = len(vector)*[0]
	x_0 = len(vector)*[1]
	k = 0

	while(dist(x_1, x_0) >= eps and k <= kmax):

		x_0 = x_1.copy()

		for i in range(len(x_1)):
			x_1[i] = vector[i]

			for j in range(i):
				x_1[i] -= matrix[i][j] * x_1[j]


			for j in range(i+1, len(x_1)):
				x_1[i] -= matrix[i][j] * x_0[j]


			x_1[i] /= matrix[i][i]

		k += 1

	print("Número de iterações: %d"%k)

	return x_0


def linear_regression(vec_x, vec_y):

    N = len(vec_x)
    sum_xy = 0
    sum_x = 0
    sum_y = 0
    sum_xx = 0

    for i in range(N):

        sum_xy += vec_x[i] * vec_y[i]
        sum_x += vec_x[i]
        sum_y += vec_y[i]
        sum_xx += vec_x[i]**2

    a1 = (sum_xy - sum_x*sum_y/N) / (sum_xx - (sum_x**2)/N)
    a0 = (sum_y - a1*sum_x)/N

    return [a0, a1]

def exponencial_adjust(vec_x, vec_y):

    new_vec_y = [ log(y) for y in vec_y ]
    coef = linear_regression(vec_x, new_vec_y)

    return [exp(coef[0]), coef[1]]

def ln_adjust(vec_x, vec_y):

    new_vec_x = [ log(x) for x in vec_x ]
    coef = linear_regression(new_vec_x, vec_y)

    return [coef[1], exp(coef[0]/coef[1])]

def hyperbolic_adjust(vec_x, vec_y):

    new_vec_x = [ (1/x) for x in vec_x ]
    coef = linear_regression(new_vec_x, vec_y)

    return coef

def potential_adjust(vec_x, vec_y):

	new_vec_y = [ log(y) for y in vec_y ]
	new_vec_x = [ log(x) for x in vec_x ]
	coef = linear_regression(new_vec_x, new_vec_y)

	return [exp(coef[0]), coef[1]]


def polinomial_adjust(vec_x, vec_y, degree):

    len_matrix = degree + 1

    #Generate Matrix
    matrix = [ [0]*(len_matrix) for i in range(len_matrix)]
    matrix[0][0] = len(vec_x)

    for j in range(1, len_matrix):
        for i in range(j+1):
            sum = 0
            for k in vec_x:
                sum += k**(i+j)

            matrix[i][j] = sum
            matrix[j][i] = sum
    #Generate Vector
    vector = []

    for i in range(len_matrix):
        sum = 0
        for k in range(len(vec_x)):
            sum += (vec_x[k]**i) * vec_y[k]
        vector.append(sum)

    return gauss_seidel(matrix, vector, 0.00001, 10000)

#Usage
vec_x = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3]
vec_y = [0.525, 0.844, 1.280, 1.863, 2.632, 3.638, 4.944, 6.626, 8.778, 11.508, 14.948]

print("Linear Regression")
print(linear_regression(vec_x, vec_y))
print("Linear Adjust")
print(polinomial_adjust(vec_x, vec_y, 1))
print("Parabolic Adjust")
print(polinomial_adjust(vec_x, vec_y, 2))
print("Exponential Adjust")
print(exponencial_adjust(vec_x, vec_y))

import matplotlib.pyplot as plt
import numpy as np

plt.scatter(vec_x, vec_y)

x = np.arange(1, 3, 0.01)
y1 = np.poly1d(linear_regression(vec_x, vec_y)[::-1])
y2 = np.poly1d(polinomial_adjust(vec_x, vec_y, 2)[::-1])
plt.plot(x, y1(x))
plt.plot(x, y2(x))

print(y2(2))

coef_exp = exponencial_adjust(vec_x, vec_y)

plt.plot(x, coef_exp[0]*np.exp(coef_exp[1]*x))


plt.show()
