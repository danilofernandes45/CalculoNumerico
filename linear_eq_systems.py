from math import sqrt

#Iterative Methods to solve a linear equations system -> O(n²)

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

def norm(x):
	sum = 0 
	for i in x:
		sum += i*i

	return sqrt(sum)

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

def gauss_jacobi(matrix, vector, eps, kmax):

	matrix = nonzero_diagonal(matrix)

	x_1 = len(vector)*[0]
	x_0 = len(vector)*[1]
	k = 0

	while(dist(x_1, x_0) >= eps and k <= kmax):

		x_0 = x_1.copy()

		for i in range(len(x_1)):
			x_1[i] = vector[i]

			for j in range(i):
				x_1[i] -= matrix[i][j] * x_0[j]


			for j in range(i+1, len(x_1)):
				x_1[i] -= matrix[i][j] * x_0[j]


			x_1[i] /= matrix[i][i]

		k += 1

	print("Numero de iterações: %d"%k)
	return x_0

#Direct Method to solve a linear equations system -> O(n³)

def schedule(matrix, vector):

 	for j in range( len(vector) - 1 ):

 		for i in range(j + 1, len(vector)):

 			m = matrix[i][j] / matrix[j][j]

 			for k in range(j, len(vector)):

 				matrix[i][k] = matrix[i][k] - m * matrix[j][k]

 			vector[i] = vector[i] - m * vector[j]

 	return matrix, vector


def solve_triangular_system(matrix, vector):

 	size = len(vector)
 	x = size * [0]

 	x[size - 1] = vector[size - 1] / matrix[size - 1][size - 1]

 	for k in range(size - 1, -1):

 		sum = 0
 		for k in range(k+1, size):
 			sum = sum + matrix[k][j] * x[j]

 		x[k] = ( vector[k] - sum) / matrix[k][k]

 	return x


def gauss_elimination(matrix, vector):

	matrix, vector = schedule(matrix, vector)
	return solve_triangular_system(matrix, vector)


matrix = [[3, 1, 1], [0, 2, 1], [0, 0, 1]]
vector = [1, 0, 1]

print("Gauss Elimination Method")
print(gauss_elimination(matrix, vector))
print()

print("Gauss-Jacobi Method")
print(gauss_jacobi(matrix, vector, 0.0001, 10000))
print()

print("Gauss-Seidel Method")
print(gauss_seidel(matrix, vector, 0.0001, 10000))
