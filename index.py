import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

fig1, ax1 = plt.subplots(1)
fig2, ax2 = plt.subplots(1)
fig3, ax3 = plt.subplots(1)
fig4, ax4 = plt.subplots(1)
fig1.tight_layout()

#solve the second order differential equation for the soliton,
def action_der(x, v, t, w):
    return (0.5*v**2 + potential_corr(x, w))

def potential(x,w): 
    return (1/4)*(x**2-1)**2 + w*(x-1)

def potential_derivative(x,w):
    return x*(x**2-1) + w

def potential_corr(x, E):
    x_plus = Newton_Raphson(2, E, 100)
    return potential(x, E) - potential(x_plus, E)

def potential_sec_der(x, w):
    return 3*x**2-1 

def Newton_Raphson(x, w, n):
    for i in range(n):
        x += -potential_derivative(x, w) / potential_sec_der(x, w)
    return x

def ODE(x, v, t, w):
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
        
    return np.array([t_total, xh_total, v_total])

def IntBisec(a_u, a_o, E, N):
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
N = 20

E = np.linspace(0.05, 0.16, N)

# Initial conditions for the Runge-Kutta algorithm.
t0 = 1e-15
v0 = 0
dt = 0.1
t_range = 60

S_Euclidean = []

def dE_dx(func, x, E):
    return np.sqrt(func(x,E))

def potential_num_w(E, x): 
    return (1/3)*x**3 - (1+2*E)*x + 2*E*(1+2*E)*np.log(1+x)

t_ends = []
phase_integral_array = []

for i in range(N):    
    
    x_range = 2
    
    a_u = -0.3
    a_o = Newton_Raphson(-1, E[i], 10)
    a_mid = IntBisec(a_u, a_o, E[i], 100)
    
    phi_mid = RK_22(a_mid, t0, v0, t_range, x_range, E[i])
    
    x_pot = np.linspace(-1.5,1.5,100)

    t = phi_mid[0]
    x = phi_mid[1]
    v = phi_mid[-1]

    #search algorithm
    n = 0
    for j in range(len(t)): 
        if abs(t[j] - 30) < 0.000001: 
            n=j
            print(n)
            break
        
    
    t_cut = t[:50]
    x_cut = x[:50]
    v_cut = v[:50]

    
    
    dB_dr = action_der(x_cut, v_cut, t_cut, E[i])
    
    ax1.plot(t_cut, x_cut, label="$\epsilon/\lambda\eta^4$ = " + str(np.round(E[i], 3)))
    ax2.plot(t_cut, dB_dr)
    
    S_e = simpson(dB_dr, dx=0.1)
    S_Euclidean.append(S_e)


ax3.plot(E, S_Euclidean, label="numerical decay amplitude")
#ax3.plot(E, 1/E, label="analytic approx.")


"""___________________________________Styles_______________________________"""
s=11
ax1.set_ylim(-1.1, 1)
#ax1.set_xlim(0, 30)
ax1.tick_params(axis='both', which='major', labelsize=s)
ax1.set_title('Several kink solutions with different $\epsilon$')
ax1.set_xlabel('$\eta\lambda^{1/2}$x', size=20)
ax1.set_ylabel('$\phi/\eta$', size=20)
ax1.grid()

ax2.set_title('Plots of $\sqrt{2V(\phi)}$ for different $\epsilon$')
ax2.set_xlabel('$\eta\lambda^{1/2}$x', size=20)
ax2.set_ylabel('$S_E\'$', size=20)
ax2.set_xlim(0,5)
ax2.grid()

ax3.set_title("Plot of minimum energy of the kink vs. $\epsilon$")
ax3.set_ylabel("$B$", size=20)
ax3.set_xlabel("$\epsilon/\lambda\eta^4 $", size = 20)
ax3.grid()

ax4.set_title("Phase space orbits of the kinks at different $\epsilon$")
ax4.set_ylabel("$\phi'$", size=20)
ax4.set_xlabel("$\phi$", size = 20)
ax4.grid()

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=s)
ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=s)

