#FINITE DIFFERENCES METHOD - ORDER 1
def dif_method(func, x, h):
	return (func(x + h) - func(x)) / h

#FINITE DIFFERENCES METHOD  - ORDER 2
def dif_advanced_method(func, x, h):
	return (-3 * func(x) + 4 * func(x + h) - func(x + 2*h)) / (2 * h)

def dif_centered_method(func, x, h):
	return (func(x + h) - func(x - h)) / (2 * h)

def dif_delayed_method(func, x, h):
	return (func(x - 2*h) - 4 * func(x - h) + 3 * func(x)) / (2 * h)


import math

print(dif_method(math.exp, 0, 0.1))
print(dif_advanced_method(math.exp, 0, 0.1))
print(dif_centered_method(math.exp, 0, 0.1))
print(dif_delayed_method(math.exp, 0, 0.1))