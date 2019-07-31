from math import exp

def trapeze_rule(func, a, b, n):

	h = (b-a)/n
	sum = 0

	for i in range(1, n):
		sum += 2 * func(a + i*h)

	sum += func(a) + func(b)

	return (h / 2) * sum

def simpson_rule(func, a, b, n):

	h = (b - a) / (2 * n)
	sum = 0

	for i in range(1, 2*n, 2):
		sum += 4 * func(a + i*h) + 2 * func(a + (i+1)*h)

	sum += func(a) - func(b)

	return (h / 3) * sum

def linear_func(x):
	return (2*x)

def gaussian_func(x):
	return exp(-x**2)

print("Linear - n = 1")
print(trapeze_rule(linear_func, 0, 1, 1))
print(simpson_rule(linear_func, 0, 1, 1))

print("Linear - n = 5")
print(trapeze_rule(linear_func, 0, 1, 5))
print(simpson_rule(linear_func, 0, 1, 5))

print("Gaussian - n = 5")
print(trapeze_rule(gaussian_func, 0, 1, 5))
print(simpson_rule(gaussian_func, 0, 1, 5))

print("Gaussian - n = 10")
print(trapeze_rule(gaussian_func, 0, 1, 10))
print(simpson_rule(gaussian_func, 0, 1, 10))

print("Gaussian - n = 100")
print(trapeze_rule(gaussian_func, 0, 1, 100))
print(simpson_rule(gaussian_func, 0, 1, 100))

print("Gaussian - n = 1000")
print(trapeze_rule(gaussian_func, 0, 1, 1000))
print(simpson_rule(gaussian_func, 0, 1, 1000))
