import numpy as np
import pandas as pd
import random
import datetime

def model_bn (n_pop, n_snp, n_indiv, frq, fst, seed: None) :
    """Simulate model bn.
        Parameters
        ----------
        n_pop : int
            Number of population(s).
        n_snp : int
            Number of SNPs
        n_indiv : int
            Number of individuals.
        frq : array
            frequency/proportion of that SNP in the data
        fst : array
            percentage of contribution to the total generic variation in each SNP
        seed : int
            Uses datetime if nothing is provided
            Randomizes generation.
        Returns
        -------
        pop_matrix : ndarray
            Individual population adimixture matrix.
        pop_idx : array
            Population adimixture matrix.

    """
    
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%% Load HapMap Info
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #HM_inf = pd.read_csv(u'hap_map',sep=u' ')
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %get allele freq and Fst for each SNP
    
    # % allele freq: frequency (proportion) of that SNP in the data
    # % Fst: % of contribution to the total genetic variation that each SNP has
    #frq = HM_inf[u'FRQ'].values #cell2mat(HM_inf(:,4));
    #Fst = HM_inf[u'FST'].values #cell2mat(HM_inf(:,3));

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
    gamma = np.zeros((n_snp,n_pop)); #allele freq for each population
    
    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    #Population matrix
    pop_matrix = np.zeros((n_pop,n_indiv)); #individual population admixture
    
    # %X = zeros(m,n); %genotype matrix

    # Uses current time if no seed is provided
    if (seed == None):
        seed = datetime.now()
    
    # % random seeding here...
    random.seed(seed)
    
    # %populate the allele freq matrix from BN with (p,F) from HapMap
    # % for each SNP...
    for i in range(0,n_snp):
        # each row will generate 'd' variates drawing from this distribution
        gamma[i,:] = np.random.beta(frq[i]*(1-fst[i])/fst[i], (1-frq[i])*(1-fst[i])/fst[i], size=n_pop)            
        
    # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

    # %populate the population admixture matrix
    # %this part is tailor made for 3 populations as described in
    # %Song et al. Nat. Genet. (2015)
    # % Treating the probabilities as ranges
    # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210
    pop_idx = np.zeros((n_indiv,1))
    
    for i in range(0,n_indiv):
        p = random.uniform(0, 1)
        if p < (60.0/210):
            pick = 1
            pop_idx[i]= 1
        elif p < (2*(60.0/210)):
            pick = 2
            pop_idx[i] = 2
        else:
            pick = 3
            pop_idx[i] = 3
            
        pop_matrix[pick-1,i] = 1
        # Leaving all other values at zero (indiv only assigned to one subpop)
    return (pop_matrix, pop_idx)
