# setup ipython environment
from ipywidgets import interact, fixed

# setup python environment
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
col = plt.rcParams['axes.color_cycle']

def mc(p, q, tmax = 1, x0 = 0):
    x = x0
    r = np.random.rand(tmax+1)
    for t in range(tmax+1):
        yield x
        if x == 0 and r[t] < 1 - p:
            x = 1
        elif x == 1 and r[t] < 1 - q:
            x = 0

def plot_mc(p:(0,1,0.1)=0.9,q:(0,1,0.1)=0.9):
    print('P=|%.1f %.1f|\n  |%.1f %.1f|' %(p,1-p,1-q,q))
    plt.figure(figsize=[9,4])
    plt.plot(list(mc(p,q,100)))
    plt.yticks([0,1],['$S_1$, off','$S_2$, on'])
    plt.xlabel('time')
    plt.ylim(-0.1,1.1)
    plt.show()

def mc_sol(s1, p, q, t):
    x0 = np.matrix([s1,1-s1])
    pij = np.matrix([[p,1-p],[1-q,q]])
    return np.array(x0*pij**t)[0]

def plot_mc_sol(s1:(0,1,0.1)=1,p:(0,1,0.1)=0.9,q:(0,1,0.1)=0.9, tmax=fixed(20)):
    s1,s2 = list(zip(*[mc_sol(s1,p,q,t) for t in range(tmax+1)]))

    plt.figure(figsize=[9,4])
    plt.plot(s1, label='S_1, off')
    plt.plot(s2, label='S_2, on')
    plt.xlabel('time')
    plt.ylabel('proportion of channels')
    plt.ylim(-0.01,1.01)
    plt.legend()
    plt.show()


def plot_mc_sol2(p:(0,1,0.1)=0.9,q:(0,1,0.1)=0.9, \
                 tmax=fixed(20), log_n:(0,10,1)=0):
    s1,s2 = list(zip(*[mc_sol(1,p,q,t) for t in range(tmax+1)]))

    n = int(np.exp(log_n))
    sim = [list(mc(p,q,tmax)) for _ in range(n)]
    sim_av = [np.mean(s) for s in list(zip(*sim))]

    print("n = %d" % n)

    plt.figure(figsize=[9,4])
    plt.plot(s2)
    plt.plot(sim_av)
    plt.xlabel('time')
    plt.ylabel('proportion of channels on')
    plt.ylim(-0.01,1.01)
    plt.show()

def red_sim(i, n, tmax=1):
    for t in range(tmax+1):
        yield i
        i = np.random.binomial(2*n,i/(2*n))

def plot_red_sim(log_n:(0,10,1)=7,prop_i:(0,1,0.1)=0.5,n_sim:(1,20,1)=1):
    n = int(np.exp(log_n))
    i = int(2*n*prop_i)
    print("n = %d, i0 = %d" % (n,i))

    plt.figure(figsize=[9,4])
    for _ in range(n_sim):
        plt.plot(list(red_sim(i,n,50)))
    plt.xlabel('time')
    plt.ylabel("number of copies of 'a'")
    plt.ylim(-2*n*0.01, 2*n*1.01)
    plt.show()

def mc4state(p1, p2, tmax=1):
    x = np.matrix([1,0,0,0])
    p = np.matrix([[1-p1/2-p2/2,p1/2,p2/2,0],[0,1-p2,0,p2],
                   [0,0,1-p1,p1],[1,0,0,0]])
    for t in range(tmax+1):
        yield x.tolist()[0]
        x = x*p

def plot_mc4state(p1:(0,1,0.1)=1, p2:(0,1,0.1)=1, \
                 plot_all=False):
    pts = list(zip(*mc4state(p1, p2, 30)))
    plt.figure(figsize=[9,4])
    plt.plot(pts[0], label='A00')
    if plot_all:
        plt.plot(pts[1], label='A10')
        plt.plot(pts[2], label='A01')
        plt.plot(pts[3], label='A11', color=col[5])
    plt.xlabel('time')
    plt.ylabel('proportion of enzymes in state')
    plt.ylim(-0.05,1.05)
    plt.legend()
    plt.show()
