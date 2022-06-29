# setup python environment
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def plot_norm(n =100):
    plt.figure(figsize=[9,4])
    x = np.linspace(-4,4)
    plt.plot(x, norm.pdf(x,0,1), lw=2)
    plt.hist(np.random.normal(size=n), normed=True, bins=30)
    plt.xlim(-4,4)
    plt.ylim(0,1)

    # show plot
    plt.show()

def plot_n_fish(n_fish):
    plt.figure(figsize=[9,4])
    plt.hist(n_fish, normed=True, bins=range(50))
    plt.xlabel('Number of fish consumed')
    plt.ylabel('Proportion of population')
    plt.xlim(0,50)
    plt.ylim(0,0.12)
    plt.show()

def plot_intakes(intakes):
    plt.figure(figsize=[9,4])
    plt.hist(intakes, normed=True, bins=30)
    plt.xlabel('Total intake (kcal)')
    plt.ylabel('Proportion of population')
    plt.xlim(0,4000)
    plt.ylim(0,0.0012)
    plt.show()
