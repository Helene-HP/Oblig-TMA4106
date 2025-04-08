import numpy as np
import matplotlib.pyplot as plt

def analytisk_f(x):
    return np.exp(x)

def nummerisk_f(x, h):
    return (np.exp(x+h) - np.exp(x))/h

def ny_formel(x, h):
    return (np.exp(x+h)-np.exp(x-h))/(2*h)

def nyny_formel(x, h):
    return (np.exp(x-2*h) - 8*np.exp(x-h) + 8*np.exp(x+h) - np.exp(x+2*h))/(12*h)


x = 1.5


print("Oppg 1:")
for i in np.logspace(0, -20, 20):
    diff = abs(analytisk_f(x) - nummerisk_f(x, i))
    print(f"h: {i:.0e}, forskjell: {diff:.8f}")

print("\nOppg 2:")
for i in np.logspace(0, -20, 20):
    diff = abs(analytisk_f(x) - ny_formel(x, i))
    print(f"h: {i:.0e}, forskjell: {diff:.8f}")

print("\nOppg 3:")
for i in np.logspace(0, -20, 20):
    diff = abs(analytisk_f(x) - nyny_formel(x, i))
    print(f"h: {i:.0e}, forskjell: {diff:.8f}")
