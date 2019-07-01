import math

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

	for i in range(len(matrix)):
		if(matrix[i][i] == 0):
			return False

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

def gauss_seidel(matrix, vector, eps, k_max):
	x_0 = len(vector)*[1]
	k = 0

	while(norm(x_1 - x_0) >= eps and k <= k_max):

		x_1 = len(vector)*[0]

		for i in range(len(x_1)):
			x_1[i] = vector[i]

			for j in range(i):
				x_1[i] -= matrix[i][j] * x_1[j]


			for j in range(i+1, len(x_1)):
				x_1[i] -= matrix[i][j] * x_0[j]


			x_1[i] /= matrix[i][i]

		x_0 = x_1

	print("Número de iterações: %d"%k)

	return x_0

#Assume-se que conhece-se o ambiente de tal forma que o retorno do ambiente e a probabilidade 	
#para cada transição de estado sob uma ação especifica são conhecidos. Estes são mapeados em matrizes assim como a política
def policy_evaluation(policy, returns, trans_prob, num_states, num_actions, amor_rate):
	#policy[state, action]
	#returns[state0, action, state1]
	#trans_prob[state0, action, state1]
	matrix = num_states * [ len(num_states) * [0] ]
	vector = len(num_states) * [0]

	for i in range(num_states):
		#state0
		matrix[i][i] = 1

		for j in range(num_actions):
			#action
			for k in range(num_states):
				#state1
				matrix[i][k] -= amor_rate * policy[i][j] * trans_prob[i, j, k]

				vector[i] += policy[i][j] * trans_prob[i, j, k] * returns[i, j, k]

	if( convergence_test(matrix) ):
		value_func = gauss_seidel(matrix, vector, eps = 0.0001, k_max = 100)
		return value_func

	return []

