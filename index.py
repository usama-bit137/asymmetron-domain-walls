import numpy as np
import matplotlib.pyplot as plt


fig1, ax1 = plt.subplots(1)
fig2, ax2 = plt.subplots(1)
fig1.tight_layout()
#solve the second order differential equation for the soliton,
def potential(x,w): 
    return (1/4)*(x**2-1)**2 + w*(x-1)

def potential_derivative(x,w):
    return x*(x**2-1)+w

def potential_sec_der(x):
    return 3*x**2-1

def Newton_Raphson(x, w, n):
    for i in range(n):
        x += -potential_derivative(x, w) / potential_sec_der(x)
    return x


def ODE(x, t, v, w):
    return potential_derivative(x, w)

def RK_22(x0, t0, v0, t_range, x_range, E):
    # Boxes to fill:
    xh_total = []
    t_total = []
    v_total = []

    while t0 < t_range:

        xh = x0 + dt * v0 / 2

        if abs(x0) > x_range:
            break

        vh = v0 + ODE(x0, v0, t0, E) * dt / 2
        x0 += dt * vh
        v0 += dt * ODE(xh, vh, t0 + dt / 2, E)
        t0 += dt

        # Fill the boxes:
        v_total.append(vh)
        xh_total.append(xh)
        t_total.append(t0)

    # Fewer complications:
    return np.array([t_total, xh_total, v_total])


# Fundamental values we require:
N = 2
E = 0

# Initial conditions for the Runge-Kutta algorithm.
t0 = 0
v0 = 0
dt = 0.1
x_range = 2
t_range = 100

x_approx = np.linspace(0,t_range,200)
y_approx = np.tanh((x_approx-t_range/2)/(np.sqrt(2)*(1-E)))

ax1.plot(x_approx,y_approx)

"""The mid-point of the overshoot and the undershoot:"""
#ax2.axvline(a_mid, color="g")

"""________________________________Finishing_Touches________________________"""
# Plots of derived action and radius:
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlim(0, t_range)

