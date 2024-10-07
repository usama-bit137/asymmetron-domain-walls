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
        x += - potential_derivative(x, w) / potential_sec_der(x)
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


def IntBisec(a_u, a_o, a_mid, E, N):
    for i in range(N):

        Phi_u = RK_22(a_u, t0, v0, t_range, x_range, E)
        amid = 0.5 * (a_u + a_o)

        Phi_mid = RK_22(amid, t0, v0, t_range, x_range, E)

        if abs(Phi_u[0, -1] - Phi_mid[0, -1]) < 1e-15:
            a_u = amid
        else:
            a_o = amid
    return amid


# Fundamental values we require:
N = 2
E = np.linspace(0.01, 0.08, 1)

# Initial conditions for the Runge-Kutta algorithm.
t0 = 0
v0 = 0
dt = 0.1
x_range = 2
t_range = 100

S = []
R = []
A_mid = []

for j in range(N):
    a_o = -0.1
    print('An overshoot value:' + str(a_o))

    a_u = Newton_Raphson(-1, E[j], 100)
    print('An undershoot value:' + str(a_u))
    
    ax2.axvline(a_o, color='r')
    ax2.axvline(a_u, color="b")
    
    """The mid-point of the overshoot and the undershoot:"""
    a_mid = 0.5 * (a_o + a_u)
    a_mid = IntBisec(a_u, a_o, a_mid, E[j], 100)
    A_mid.append(a_mid)
    Phi_mid = RK_22(a_mid, t0, v0, t_range, x_range, E[j])

    # Inputs for the action:
    t = Phi_mid[0, :]
    x = Phi_mid[1, :]
    v = Phi_mid[-1, :]
    
    # Removing the waste end:
    for l in np.arange(0, len(t) - 1):
        if np.round(Phi_mid[0, l]) == 40:
            n = round(l)
            t_red = t[:n]
            x_red = x[:n]
            v_red = v[:n]
            break
    
    
        ax1.plot(t,x)
        ax2.plot(Phi_mid[1,:], -potential(Phi_mid[1,:], E[j]))
        ax2.axvline(a_o, color='r')
"""________________________________Finishing_Touches________________________"""
# Plots of derived action and radius:
ax1.set_ylim(-1.5, 1.5)
ax1.set_xlim(5, 100)

