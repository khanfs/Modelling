# setup ipython environment
from ipywidgets import interact, fixed

# setup python environment
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
plt.style.use('ggplot')
col = plt.rcParams['axes.color_cycle']

def plot_exp(double_time:(5,60,5)=20, ic=fixed(100), tmax=fixed(60)):
    K = np.log(2)/double_time

    t = np.linspace(0, tmax, 30)

    plt.figure(figsize=[9,4])
    plt.plot(t, ic*np.exp(K*t), lw=2)
    plt.plot([double_time,double_time,0],[0,ic*2,ic*2], lw=2)
    plt.xlabel('time')
    plt.ylabel('population density')
    plt.xlim(0,tmax)
    plt.show()

    print("K = %.3f" % K)

def plot_lf(r:(0,4,0.1), x0=fixed(0.05), tmax=fixed(10)):
    t = np.linspace(0, tmax, 50)
    plt.figure(figsize=[9,4])
    plt.plot(t, x0/(x0+(1-x0)*np.exp(-r*t)), lw=2)
    plt.axhline(1, ls='--', lw=1, c='k')
    plt.xlabel('time')
    plt.ylabel('population density')
    plt.xlim(0,tmax)
    plt.ylim(0,1.05)
    plt.show()

def lv(xy, t, r):
    x,y = xy
    return [x*(1 - y), y*(x - r)]

def plot_lv(r:(0.1,5,0.1)=0.3, x0:(0,10,1)=10, y0=fixed(1)):
    t = np.linspace(0.0, 20.0, 200)
    x,y = list(zip( *integrate.odeint(lv, [x0,y0], t, (r,)) ))

    plt.figure(figsize=[9,4])
    plt.plot(t, x, lw=2, label="prey")
    plt.plot(t, y, lw=2, label="predator")
    plt.xlabel('time')
    plt.ylabel('population density')
    plt.ylim(0, 14)
    plt.legend()
    plt.show()

def plot_lv_pp(r=fixed(1.5), x0:(0.4,4,0.1)=2, y0=fixed(1)):
    t = np.linspace(0.0, 5.0, 30)
    x,y = list(zip( *integrate.odeint(lv, [x0,y0], t, (r,)) ))
    x = np.array(x); y = np.array(y)

    plt.figure(figsize=[9,6])

    plt.quiver(x[:-1], y[:-1], \
               (x[1:]-x[:-1])*0.98, (y[1:]-y[:-1])*0.98, \
               scale_units='xy', angles='xy', scale=1, \
               width=0.005, headwidth=3, headlength=3, \
               headaxislength=3, color=col[0])

    plt.scatter([r, x0], [1, y0], s=60, marker='o', c='k')
    plt.xlabel('prey population')
    plt.ylabel('predator population')
    plt.xlim(0, 4.2)
    plt.ylim(0,3.4)
    plt.show()

def plot_lv_stream(r=1.5):
    gx, gy = np.meshgrid(np.linspace(0, 4.21, 11), np.linspace(0, 3.41, 11))
    dx, dy = lv([gx, gy], 0, r)

    plt.figure(figsize=[9,6])
    plt.streamplot(gx, gy, dx, dy, color=col[0])
    plt.scatter(r, 1, s=60, marker='o', c='k')

    plt.xlabel('prey population')
    plt.ylabel('predator population')
    plt.xlim(0, 4.2)
    plt.ylim(0,3.4)
    plt.show()

# 3D example, Lorenz
from mpl_toolkits.mplot3d import Axes3D

def lorenz_attractor(xyz, t, S = 10.0, R = 28.0, B = 8.0/3.0):
    x, y, z = xyz
    return [S*(y - x), -x*z + R*x - y, x*y - B*z]

la_pts = list(zip(* integrate.odeint(lorenz_attractor, [1, 1, 1],
                                     np.linspace(0,100,1e5))[100:] ))

def plot_la(angle:(0,360,10)=0):
    # get onto attractor

    fig = plt.figure(figsize=[9,4])
    ax = fig.add_subplot(111, projection='3d', axisbg="white")
    ax.view_init(elev=10., azim=angle)
    ax.plot(*la_pts, lw=0.1, label="prey")
    ax.axis('off')

    plt.show()

def plot_ss(beta:(-4,4,0.5)=-1, gamma:(-4,4,0.5)=1):
    def f(xy, beta, gamma):
        x,y=xy
        a11 = 1; a12 = 2
        a22 = beta - a11
        a21 = (a11*a22 - gamma)/a12
        return [a11*x + a12*y, a21*x + a22*y]

    l1 = 0.5*(beta+np.sqrt(beta-4*gamma+0j))
    l2 = 0.5*(beta-np.sqrt(beta-4*gamma+0j))
    print("eigenvalues: %.3f + %.3fi, %.3f + %.3fi" \
            % (l1.real, l1.imag, l2.real, l2.imag))

    gx, gy = np.meshgrid(np.linspace(-1, 1, 11), \
                         np.linspace(-1, 1, 11))
    dx, dy = f([gx, gy], beta, gamma)

    fig = plt.figure(figsize=[9,6])
    ax = fig.add_subplot(111)
    plt.streamplot(gx, gy, dx, dy)
    plt.scatter(0, 0, s=60, marker='o', c='k')
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    ax.axis('off')
    plt.show()

def sir(si, t, beta, gamma, nu):
    s,i = si
    return [-beta*s*i + gamma*(1 - s - i), \
            beta*i*s - nu*i]

def plot_sir(sigma:(0,2,0.1)=1, gamma:(0,1,0.1)=0.1):
    nu = 1; beta = sigma*nu

    t = np.linspace(0.0, 50.0, 100)
    s,i = list(zip( *integrate.odeint(sir, [0.9,0.1], t,\
                                      (beta, gamma, nu)) ))

    plt.figure(figsize=[9,4])
    plt.plot(t, s, label="susceptible")
    plt.plot(t, i, label="infected")
    plt.plot(t, [1-si-ii for si,ii in zip(s,i)], label="recovered")
    plt.xlabel('time')
    plt.ylabel('population density')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()

def plot_sir_pp(sigma:(0.25,2,0.25)=0.5, gamma:(0,1,0.5)=0.5):
    nu = 1; beta = sigma*nu

    t = np.linspace(0.0, 50.0, 100)
    s,i = list(zip( *integrate.odeint(sir, [0.9,0.1], t,\
                                      (beta, gamma, nu)) ))

    gx, gy = np.meshgrid(np.linspace(0, 1, 11), \
                         np.linspace(0, 1, 11))
    dx, dy = sir([gx, gy], 0, beta, gamma, nu)

    plt.figure(figsize=[9,6])
    plt.streamplot(gx, gy, dx, dy)

    plt.scatter([1, nu/beta], [0, gamma*(beta-nu)/(beta*(nu+gamma))], \
                s=60, marker='o', c='k')

    plt.xlabel('susceptible')
    plt.ylabel('infected')
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.show()
