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

			x_1[i] /= matrix[i][i]

		k += 1

	print("Número de iterações: %d"%k)

	return x_0				x_1[i] -= matrix[i][j] * x_0[j]


			x_1[i] /= matrix[i][i]

		k += 1

	print("Número de iterações: %d"%k)

	return x_0


def exponencial_adjust(vec_x, vec_y):

    new_vec_y = [ ln(y) for y in vec_y ]
    coef = linear_regretion(vec_x, new_vec_y)

    return [exp(coef[0]), coef[1]]


def linear_regretion(vec_x, vec_y):

    N = len(vec_x)
    sum_xy = 0
    sum_x = 0
    sum_y = 0
    sum_xx = 0

    for i in range(N):

        sum_xy = vec_x[i] * vec_y[i]
        sum_x = vec_x[i]
        sum_y = vec_y[i]
        sum_xx = vec_x[i]**2

    a1 = (sum_xy - sum_x*sum_y/N) / (sum_xx - (sum_x**2)/N)
    a0 = (sum_y - a1*sum_x)/N

    return [a0, a1]

def polinomial_adjust(vec_x, vec_y, degree):

    len_matrix = degree + 1

    #Generate Matrix
    matrix = [ [0]*(len_matrix) for i in range(len_matrix)]
    matrix[0][0] = len(vec_x)

    for j in range(1, len_matrix):
        for i in range(j+1):
            sum = 0
            for k in vec_x:
                sum += vec_x**(i+1)

            matrix[i][j] = matrix[j][i] = sum

    #Generate Vector
    vector = []

    for i in range(len_matrix):
        sum = 0
        for k in range(len(vec_x)):
            sum += (vec_x[k]**i) * vec_y[k]
        vector.append(sum)

    return gauss_seidel(matrix, vector, 0.00001, 10000)
