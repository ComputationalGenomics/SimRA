import numpy as np
import scipy as sp
import pandas as pd
import sys
import math
import random

def simulate(v_noise, k, Pn, Vn, N, prop1, div=2, seed=0):
    """
    Returns a matrix and a dictionary of phenotypes associated to different variables.

    Parameters
    ----------
    v_noise: float
        Noise tolerance
    k: list
        The number of combinating variables contributing to each phenotype
    Pn: int 
        The number of phenotypes
    Vn: int
        The number of variables
    N: int 
        The number of samples
    prop1: float
        Proportion of one's in each vector zero matrix
    div: int, optional
    seed: int, optional
        Sets the seed value for the random seed() method
    """
    np.random.seed(seed)
    random.seed(seed)

    numzero = int(prop1*N)

    varlist = []
    # generate random variables
    for c in range(Vn):
        arr = np.array([0]*numzero + [1]*(N-numzero))
        np.random.shuffle(arr)
        varlist.append(arr)

    varmat = np.vstack(varlist)
    # print(varmat)

    propK = np.array(list(k[0]*np.ones(div)) + list(k[1]*np.ones(Pn-div)))
    # print(propK)
    propK = propK.astype(int)

    # obtain number of interacting variables
    #propK = random.sample(range(2, k), Pn)

    Kmap = []
    phenodict = {}

    for p in range(Pn):
        randidx = np.random.choice(varmat.shape[0], propK[p], replace=False)
        Kmap.append(randidx)
        randmat = varmat[randidx, :]

        #pheno1 = randmat.all(axis=0)
        pheno = np.logical_and.reduce(randmat, axis=0)
        pheno = pheno.astype(int)
        if sum(pheno) == 0:
            print(randidx)

        # introduce random noise for 1-v_bin
        flipidx = random.sample(range(0, N), int(v_noise*N))
        for i in flipidx:
            pheno[i] = 1 - pheno[i]

        phenodict[tuple(randidx)] = pheno
    return phenodict, varmat