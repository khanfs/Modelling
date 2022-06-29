# setup ipython environment
from ipywidgets import interact, fixed

# setup python environment
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.style.use('ggplot')
col = plt.rcParams['axes.color_cycle']


def hh(x, t, I=0):
    v, n, m, h = x
    c = 1; gna = 120; gk = 36; gl = 0.3
    vna = 115; vk = -12; vl = 10.5989

    an = lambda v: 0.01*(-v+10)/(np.exp((-v+10)/10)-1)
    bn = lambda v: 0.125*np.exp(-v/80)
    am = lambda v: 0.1*(-v+25)/(np.exp((-v+25)/10)-1)
    bm = lambda v: 4*np.exp(-v/18)
    ah = lambda v: 0.07*np.exp(-v/20)
    bh = lambda v: 1/(np.exp((-v+30)/10)+1)

    dv = -(gna*m**3*h*(v - vna) + gk*n**4*(v - vk) + gl*(v - vl) - I)/c
    dn = an(v)*(1 - n) - bn(v)*n
    dm = am(v)*(1 - m) - bm(v)*m
    dh = ah(v)*(1 - h) - bh(v)*h

    return [dv, dn, dm, dh]

def plot_hh(Iapp:(0,30,2.5)=0):

    t = np.linspace(0,100,1000)
    x = list(zip(* integrate.odeint(hh, [0,0,0,0], t, (Iapp,)) ))

    plt.figure(figsize=[9,6])
    plt.plot(t,np.array(x[0])-60)
    plt.xlabel('t')
    plt.ylabel('V')
    plt.ylim(-80,40)
    plt.show()

def fitzn(x, t, I, a, b, c):
    v, w = x
    return [v*(a - v)*(v - 1) - w + I, b*v - c*w]

def plot_fitzn(v0:(0,1,0.1)=0,Iapp:(0,1,0.1)=0, a:(0,1,0.1)=0.25,
               b:(0,0.5,0.002)=0.002, c:(0,0.5,0.002)=0.002):

    t = np.linspace(0,1e3,200)
    x = list(zip(* integrate.odeint(fitzn, [v0,0], t, (Iapp,a,b,c)) ))

    plt.figure(figsize=[9,6])
    plt.plot(t,x[0], label='v')
    plt.plot(t,x[1], label='w')
    plt.xlabel('t')
    plt.ylabel('V')
    plt.ylim(-0.5,2)
    plt.legend()
    plt.show()

def plot_fitzn_pp(v0:(0,1,0.1)=0,Iapp:(0,1,0.1)=0, a:(0,1,0.1)=0.25,
               b:(0,0.5,0.002)=0.002, c:(0,0.5,0.002)=0.002):
    gv, gw = np.meshgrid(np.linspace(-0.5, 1.5, 15), \
                         np.linspace(-0.2, 1.0, 15))
    dv, dw = fitzn([gv, gw], 0, Iapp, a, b, c)

    nullv = lambda v: v*(a - v)*(v - 1) + Iapp
    nullw = lambda w: b/c*v

    v = np.linspace(-4,4,100)
    w1 = nullv(v)
    w2 = nullw(v)

    t = np.linspace(0,1e3,200)
    x = list(zip(* integrate.odeint(fitzn, [v0,0], t, (Iapp,a,b,c)) ))

    plt.figure(figsize=[9,6])
    plt.streamplot(gv, gw, dv, dw, color='grey')
    plt.plot(v,w1,color='k')
    plt.plot(v,w2,color='k')
    plt.plot(x[0],x[1],color=col[0], lw=3)
    plt.scatter(v0,0, s=60, marker='o', color=col[0], zorder=10)
    plt.xlabel('v')
    plt.xlim(-0.5,1.5)
    plt.ylabel('w')
    plt.ylim(-0.2,1.0)
    plt.show()

def switch(p,t,kact_s,kdecay,k,n):
    return kact_s + p**n/(k**n + p**n) - kdecay*p

def plot_switch(p0:(0,10,1)=0, \
                kact_s:(0,0.2,0.01)=0.15, kdecay:(0,0.3,0.05)=0.1, \
                k:(0.5,5,0.4)=4, n:(1,5,1)=3):
    t = np.linspace(0,100,100)
    p = integrate.odeint(switch,p0,t,(kact_s,kdecay,k,n))

    plt.figure(figsize=[9,6])
    plt.plot(t,p)
    plt.xlabel('t')
    plt.ylabel('p')
    plt.ylim(0,12)
    plt.show()


def plot_switch_eqns(kact_s:(0,1,0.05)=0.15, kdecay:(0,0.3,0.05)=0.1, \
                 k:(0.5,5,0.4)=4, n:(1,5,1)=3):
    nc0 = lambda p: p**n/(k**n + p**n)
    nc1 = lambda p: kdecay*p - kact_s

    p_pts = np.linspace(0,12,30)

    plt.figure(figsize=[9,6])
    plt.plot(p_pts,nc0(p_pts))
    plt.plot(p_pts,nc1(p_pts))
    plt.xlabel('p')
    plt.ylim(-0.2,1.1)
    plt.show()

def plot_switch_ss(kdecay:(0,0.2,0.01)=0.2, k:(0.5,5,0.6)=4): # n=2
    coeff = lambda kact_s: [-kdecay, kact_s + 1, \
                            -kdecay*k**2, kact_s*k**2]

    k_pts = np.linspace(0,0.3,500)
    ss = list(zip(* [np.roots(coeff(ka)) for ka in k_pts]))
    ss_mask = [np.imag(ssi)==0 for ssi in ss]
    ss_re = [np.array(np.real(ssi))[mask] for ssi, mask in zip(ss,ss_mask)]

    plt.figure(figsize=[9,6])
    plt.plot(k_pts[ss_mask[0]], ss_re[0], color=col[0])
    plt.plot(k_pts[ss_mask[1]], ss_re[1], color=col[0], ls='--')
    plt.plot(k_pts[ss_mask[2]], ss_re[2], color=col[0])
    plt.xlabel('k_act s')
    plt.xlim(0,0.3)
    plt.ylabel('p')
    plt.ylim(0,12)
    plt.show()
