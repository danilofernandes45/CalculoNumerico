from math import sqrt

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

#O(n²)
def mne_vandermonde(vec_x, vec_y):

    matrix = []
    for i in range(len(vec_x)):
        row = [vec_x[i]**j for j in range(len(vec_x))]
        matrix.append(row)

    return gauss_seidel(matrix, vec_y, 0.000001, 10000)

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


vec_x = [3, 7, 14, 21, 28]
vec_y1 = [14.4025, 27.425, 34.49, 34.375, 36.165]
vec_y2 = [21.845, 31.28, 35.16, 38.1, 39.98]
vec_y3 = [19.0045, 26.745, 31.3, 34.18, 35.665]

print("Result from Vandermond Method")
print(vandermond_method(vec_x, vec_y1))

print("Result from Newton Method (Operators)")
print(newton_method(vec_x, vec_y1))

print("Result from Vandermond Method")
print(vandermond_method(vec_x, vec_y2))

print("Result from Newton Method (Operators)")
print(newton_method(vec_x, vec_y2))

print("Result from Vandermond Method")
print(vandermond_method(vec_x, vec_y3))

print("Result from Newton Method (Operators)")
print(newton_method(vec_x, vec_y3))
