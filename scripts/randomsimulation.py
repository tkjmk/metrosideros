#!/usr/bin/python3 -W ignore
__author__ = 'Tane Kafle (tk916@ic.ac.uk)'
__version__ = '0.0.2'
__date__ = 'May 2018'

"""
Changelog:
----------
v0.0.2 - May 2025
    - Refactored to include main() function and argparse handling

v0.0.1 - May 2018
    - Initial version
"""

import numpy as np
import matplotlib.pyplot as plt

def getrandom(samplesize=100, prob1=1.0/2, prob2=1.0/2):
    rdm = sum(np.random.choice([0, 1], size=(samplesize,), p=[prob1, prob2]))
    return rdm
# The length returned is the number of positives assigned.

def main():
    g = [None] * 1000 #Change this for list of randomisation tend
    for i in range(0,len(g)):
        g[i] = getrandom(670, 0.5,0.5)
        
    sum(g)/len(g)


    g = np.array(g)
    m = np.mean(g)
    s = np.std(g)
    p = np.arange(1)


    fig, ax = plt.subplots()
    ax.bar(p, m, yerr=s, align='center', alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylabel('Positives')
    ax.set_xticks(p)
    ax.set_xticklabels("Random")
    ax.set_title('Model')
    ax.yaxis.grid(True)

if __name__ == "__main__":
    main()
