# setup ipython environment
from ipywidgets import interact, fixed

# setup python environment
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def plot_expm(r:(0,4,0.1)=1.1, n0=fixed(100)):
    nt = [r**t * n0 for t in range(10)]

    plt.figure(figsize=[9,4])
    plt.plot(nt, lw=2)
    plt.xlabel('generation number')
    plt.ylabel('population size')
    plt.ylim(0, max(nt)*1.1)
    plt.show()

def logistic_map(xt, r):
     return r*xt*(1-xt)

def gen_lm(xt, r, tmax=1):
    for t in range(tmax+1):
        yield xt
        xt = logistic_map(xt, r)

def plot_lm(r:(0,2,0.1)=1.1):
    plt.figure(figsize=[9,4])
    for ic in np.linspace(0.1, 0.9, 5):
        plt.plot(list(gen_lm(ic, r, 20)), lw=2)
    plt.xlabel('generation number')
    plt.ylabel('population size')
    plt.ylim(0,1)
    plt.show()

def plot_lm_web(r:(0,2,0.05)=1.1, ic:(0,1,0.01)=0.5, n_iter:(1,100)=1):
    plt.figure(figsize=[8,7])
    plt.plot([0,1], [0,1], 'k') # x=y line

    # logistic map
    x_pts = np.linspace(0,1,30)
    plt.plot(x_pts, logistic_map(x_pts, r), lw=2)

    # plot 'web'
    lm_pts = [v for v in gen_lm(ic, r, n_iter) for _ in (0, 1)]
    plt.plot(lm_pts[:-2], [0]+lm_pts[2:-1], lw=1)

    plt.xlabel('x_t')
    plt.ylabel('x_{t+1}')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.show()

def plot_lm_bif1():
    plt.figure(figsize=[9,4])
    c = plt.rcParams['axes.color_cycle']

    r_pts1 = np.linspace(0.001, 1, 30)
    r_pts2 = np.linspace(1,     3, 30)
    plt.plot(r_pts1, [0]*len(r_pts1), ls='-',  lw=2, color=c[2])
    plt.plot(r_pts2, [0]*len(r_pts2), ls='--', lw=2, color=c[2])
    plt.plot(r_pts1, 1-1/r_pts1, ls='--', lw=2, color=c[5])
    plt.plot(r_pts2, 1-1/r_pts2, ls='-',  lw=2, color=c[5])

    plt.xlabel('$r$')
    plt.ylabel('$x*$')
    plt.ylim(-1,1)
    plt.show()

def plot_lm_bif2():
    def lm_ss(x0, r):
        gen = gen_lm(0.1, r, 1000)
        [next(gen) for _ in range(800)] #ignore first 800 pts
        return list(set([round(pt,3) for pt in gen]))

    def lm_unique_points_to_plot(r_pts):
        for r in r_pts:
            for pt in lm_ss(0.1, r):
                yield (r, pt)

    r_pts  = np.linspace(2.8,4,1000)
    lm_pts = list(zip(*lm_unique_points_to_plot(r_pts)))

    plt.figure(figsize=[9,5])
    plt.plot(lm_pts[0], lm_pts[1], ',')
    plt.xlabel('$r$')
    plt.xlim(2.8,4)
    plt.ylabel('$x*$')
    plt.show()

# plant model
def gen_pm(a, b, g, s, p0, p1, tmax=2):
    p = [p0, p1]
    yield p0; yield p1
    for t in range(tmax-1):
        pn = a*s*g*p[1] + b*s**2*(1-a)*g*p[0]
        yield pn
        p = [p[1], pn]

def plot_pm(a:(0,1,0.05)=1, b:(0,1,0.05)=1, g:(0,2,0.1)=1, s:(0,1,0.05)=1):
    plt.figure(figsize=[8,7])
    pts = list(gen_pm(a, b, g, s, 0, 10, 20))
    plt.plot(pts, lw=2)
    plt.xlabel('generation number')
    plt.ylabel('population size')
    plt.ylim(0, max(pts)*1.1)
    plt.show()
