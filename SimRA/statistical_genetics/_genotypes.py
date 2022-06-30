import numpy as np
import pandas as pd
import random
import datetime

def model_bn(hap_map, d, m, n, seed=None):
    """Simulate model bn.
        Parameters
        ----------
        hap_map : a dataframe
            Training data.
        d : integer
            Number of population(s).
        m : integer
            Number of SNPs
        n : integer
            Number of individuals.
        seed : number
            Uses datetime if nothing is provided
            Randomizes generation.
    """
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%% Load HapMap Info
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HM_inf = pd.read_csv(u'hap_map',sep=u' ')
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %get allele freq and Fst for each SNP
    
    # % allele freq: frequency (proportion) of that SNP in the data
    # % Fst: % of contribution to the total genetic variation that each SNP has
    frq = HM_inf[u'FRQ'].values #cell2mat(HM_inf(:,4));
    Fst = HM_inf[u'FST'].values #cell2mat(HM_inf(:,3));

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # % REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # % HapMap3 Data ~1000 individuals and ~1000000 SNPs
    # m = int(1e4) #number of SNPs
    # n = int(1e3) #number of individuals
    
    # % each row of Gamma will be populated with 3 i.i.d. draws from BN model
    # % thanks to the nature of the HapMap data (sampled from 3 populations)
    #
    # % allele freq of population: allele freq of each SNP described by that
    # % population
    G = np.zeros((m,d)); #allele freq for each population
    
    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    S = np.zeros((d,n)); #individual population admixture
    
    # %X = zeros(m,n); %genotype matrix

    # Uses current time if no seed is provided
    if (seed == None):
        seed = datetime.now()
    
    # % random seeding here...
    random.seed(seed)
    
    # %populate the allele freq matrix from BN with (p,F) from HapMap
    # % for each SNP...
    for i in xrange(0,m):
        # each row will generate 'd' variates drawing from this distribution
        G[i,:] = np.random.beta(frq[i]*(1-Fst[i])/Fst[i], (1-frq[i])*(1-Fst[i])/Fst[i], size=d)            
        
    # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

    # %populate the population admixture matrix
    # %this part is tailor made for 3 populations as described in
    # %Song et al. Nat. Genet. (2015)
    # % Treating the probabilities as ranges
    # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210
    popidx = np.zeros((n,1));
    
    for i in xrange(0,n):
        p = random.uniform(0, 1);
        if p < (60.0/210):
            pick = 1;
            popidx[i]= 1;
        elif p < (2*(60.0/210)):
            pick = 2;
            popidx[i] = 2;
        else:
            pick = 3;
            popidx[i] = 3;
            
        S[pick-1,i] = 1;
        # Leaving all other values at zero (indiv only assigned to one subpop)
