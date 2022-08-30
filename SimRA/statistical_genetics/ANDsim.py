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


def generate(Vn, Pn, varmat, phenodict):
    """
    Names the variables generated in the simulate function and converts the dictionary of phenotypes into a dataframe.
    Writes the dataframes to an excel file.

    Parameters
    ----------
    Pn: int 
        The number of phenotypes
    Vn: int
        The number of variables
    varmat: array
        Matrix of random variables
    phenodict: dictionary
        Contributing variables mapped to each phenotype
    """
    # name variables
    varlist = ['var'+str(i) for i in range(0, Vn)]
    phenolist = ['pheno'+str(i) for i in range(0, Pn)]

    Vdf = pd.DataFrame(varmat.T, columns=varlist)

    Pdf = pd.DataFrame.from_dict(phenodict)
    Pdf.columns = phenolist

    SimDF = pd.concat([Pdf, Vdf], axis=1)

    # convert phenodict into dataframe
    idx = 0
    mapping = {}
    for key in phenodict.keys():
        varkey = [varlist[j] for j in key]
        mapping[phenolist[idx]] = [varkey]
        idx += 1

    mapDF = pd.DataFrame.from_dict(mapping).T
    mapDF.columns = ['Contributing Variables']

    dirpath = sys.path[0]
    outfname = dirpath+"/Sim_eps"+str(v)+"_k"+str(k)+"_And.xlsx"

    with pd.ExcelWriter(outfname) as writer:
        SimDF.to_excel(writer, sheet_name='Data')
        mapDF.to_excel(writer, sheet_name='Mapping')


def run_simulation(v_noise, k, Pn, Vn, N, prop1, seed):
    """
    Allows the use of the simulate and generate function at the same time

    Parameters
    ----------
    Same as the simulate() function
    """
    phenodict, varmat = simulate(v_noise, k, Pn, Vn, N, prop1, seed)
    generate(Vn, Pn, varmat, phenodict)
    return ("End of simulation")