import numpy as np
from scipy.optimize import root

def adam_moulton_fourth_order(f, y0, t):
    num_steps = len(t)
    y = np.zeros((num_steps, len(y0)))
    y[0] = y0

    def predictor_corrector(y_i, t_i, h):
        y_ip1=y_i+h/24*(55*f(t_i,y[i])-59*f(t_i,y[i-1])+37*f(t_i,y[i-2]) - 9*f(t_i,y[i-3]))
        #for _ in range(2):
        y_ip1= y_i + h/24 * (9 * f(t_i + h, y_ip1) + 19 * f(t_i, y_i) - 5 * f(t_i - h, y[i - 1]) + f(t_i - 2 * h, y[i - 2]))
        return y_ip1

    for i in range(0, num_steps - 1):
        y_i = y[i]
        t_i = t[i]
        h = t[i + 1] - t_i

        if i < 3:
            # Use a lower-order method (e.g., Runge-Kutta) for the first few steps
            k1 = h * f(t_i, y_i)
            k2 = h * f(t_i + h / 2, y_i + k1 / 2)
            k3 = h * f(t_i + h / 2, y_i + k2 / 2)
            k4 = h * f(t_i + h, y_i + k3)
            y[i + 1] = y_i + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        else:
            y[i + 1] = predictor_corrector(y_i, t_i, h)

    return y

# Example usage:
# Define your ODE function f(t, y)
def f(t, y):
    return y*(1.+0.5/t)
def sol(t):
    return np.sqrt(t)*np.exp(t)
# Define the time steps and initial condition
t = np.linspace(1., 2., 1000)
y0 = np.array([np.exp(1.0)])

# Solve the ODE using Adam-Moulton (4th order)
result = adam_moulton_fourth_order(f, y0, t)
for i in range(1000):
    print((result[i] - sol(t[i]))/sol(t[i]))
