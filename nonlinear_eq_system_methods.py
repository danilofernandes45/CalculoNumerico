#Methods to resolution of Non-linear Equations System
from math import sqrt

def dist(x1, x2):
	sum = 0
	for i in range(len(x1)):
		sum += (x1[i] - x2[i])**2

	return sqrt(sum)

def norm(x):
	sum = 0
	for i in x:
		sum += i ** 2
	return sqrt(sum)

def sum(x1, x2):
	y = len(x1) * [0]
	for i in range(len(x1)):
		y[i] = x1[i] + x2[i]

	return y

def neg(x):
	for i in range(len(x)):
		x[i] = -1 * x[i]

	return x

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

#Method to solve linear equation system
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

def function(x):
	y1 = x[0]**2 + x[1]
	y2 = x[0] * x[1]

	return [ y1, y2 ]

def jacobian(x):
	j11 = 2*x[0]
	j12 = 1
	j21 = x[1]
	j22 = x[0]

	return [[j11, j12], [j21, j22]]

def newton_rapson(function, jacobian, x, eps, kmax):

	k = 0

	delta_x = len(x)*[1]

	while(norm(delta_x) >= eps and k < kmax):
		
		delta_x = gauss_seidel(jacobian(x), neg(function(x)), eps, kmax)
		x = sum(x, delta_x)

		k += 1

	print("Número de iterações: %d"%k)
	return x

def modified_newton_rapson(function, jacobian, x, eps, kmax):

	k = 0

	delta_x = 1

	jacobian = jacobian(x)

	while(delta_x >= eps and k < kmax):
		
		delta_x = gauss_seidel(jacobian, neg(function(x)), eps, kmax)
		x = x + delta_x

		k += 1

	print("Número de iterações: %d"%k)
	return x


#Usage

matrix = [[3, 1], [1, 2]]
vector = [1, 0]

print(gauss_seidel(matrix, vector, 0.0001, 1000))

print(nonzero_diagonal([[0, 1, 2], [1, 0, 2], [0, 0, 1]]))

print(newton_rapson(function, jacobian, [0, 0], 0.0001, 1000))