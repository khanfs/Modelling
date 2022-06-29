# setup ipython environment
from ipywidgets import interact, fixed

# setup python environment
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.style.use('ggplot')
col = plt.rcParams['axes.color_cycle']

def diffusion(log_D:(-20,0,1)=-5):
    D = 10**log_D # units of D are cm^2/sec
    dist = [10**log_d for log_d in range(-6,3)]
    t = [d**2/D for d in dist]

    plt.figure(figsize=[9,6])
    plt.plot(dist,t)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(0.001,60*60)
    # units are cm
    plt.xticks([1e-4,1e-3,1e-1,1e0,1e2], \
               ['1 um','10 um','1mm','1cm','1m'])
    plt.yticks([0.001,1,60,60*60,60*60*24,60*60*24*365,60*60*24*365*100], \
               ['0.001 sec','1 sec','1 min','1 hour','1 day','1 year','100 years'])
    plt.show()

def diff2d(L:(40,500,10)=280,a:(5,40,5)=20,D:(1,40,0.5)=22.5):
    tau = L**2/(2*D)*np.log(L/a)
    print("tau = %d min = %d hours" % (tau, round(tau/60)))

def fisher(uv, t, c):
    u,v = uv
    return [v, -c*v-u*(1-u)]

def plot_fisher(c:(0,4,0.1)=1):
    gx, gy = np.meshgrid(np.linspace(-0.5, 1.5, 11), \
                         np.linspace(-1, 1, 11))
    dx, dy = fisher([gx, gy], 0, c)

    plt.figure(figsize=[9,6])
    plt.streamplot(gx, gy, dx, dy, color=col[0])
    plt.xlabel('U')
    plt.ylabel('V')
    plt.xlim(-0.5,1.5)
    plt.ylim(-1,1)
    plt.show()

def plot_fisher2(c:(0,4,0.1)=1):
    gx, gy = np.meshgrid(np.linspace(0, 1.1, 11), \
                         np.linspace(-0.5, 0.1, 21))
    dx, dy = fisher([gx, gy], 0, c)

    t = np.linspace(0, 100, 500)
    u,v = list(zip(* integrate.odeint(fisher, [0.9999,0], t, (c,)) ))

    plt.figure(figsize=[9,6])
    plt.streamplot(gx, gy, dx, dy, color=col[0])
    plt.scatter([0,1],[0,0], s=60, marker='o', c='k')
    plt.plot(u,v, c='k')
    plt.xlabel('U')
    plt.ylabel('V')
    plt.xlim(-0.01,1.01)
    plt.ylim(-0.5,0.01)
    plt.show()

def cell(rule_no, x, tmax=1, wrap=True):
    keys = [(1,1,1),(1,1,0),(1,0,1),(1,0,0),
            (0,1,1),(0,1,0),(0,0,1),(0,0,0)]
    vals = [int(v) for v in list('{:08b}'.format(rule_no))]
    rules = dict(zip(keys,vals))
    print(rules)

    for t in range(tmax+1):
        yield x
        x = [x[-1]] + x + [x[1]]
        x = [rules[xt] for xt in zip(x, x[1:], x[2:])]

def plot_cell(rule_no:(0,255,1)=0):
    np.random.seed(10)
    x0 = list(np.random.randint(0,2,100))
    gol = list(cell(rule_no,x0,100))

    plt.figure(figsize=[8,8])
    plt.imshow(gol,interpolation='nearest', cmap=cm.binary)
    plt.axis('off')
    plt.show()

def dcell_cycle(x,t,l=0.1):
    return [l*(2*x[-1]-x[0])] + list(l*(x[:-1] - x[1:]))

def plot_dcell_cycle(l:(0,0.3,0.05)=0,k:(1,10,1)=3):
    t = np.linspace(0,100,100)
    x0 = [1] + [0]*(k-1)
    x = integrate.odeint(dcell_cycle,x0,t,(l,))

    plt.figure(figsize=[9,6])
    plt.plot(t,x)
    plt.xlabel('t')
    plt.ylabel('N_i')
    plt.ylim(0,5)
    plt.show()
