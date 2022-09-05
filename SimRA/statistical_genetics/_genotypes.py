import numpy as np
import pandas as pd
import random
import datetime

def model_bn (n_pop, n_snp, n_indiv, frq, fst, seed: None):
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
            frequency/proportion of a specific SNP in the data
        fst : array
            percentage of contribution to the total generic variation from each SNP
        seed : int
            Uses datetime if nothing is provided
            Randomizes generation.
        Returns
        -------
        pop_matrix : 2d data array
            Individual population adimixture matrix.
        pop_idx : array
            Population adimixture matrix.
        genetic_matrix : 2d data array
            Genetic data matrix. 

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
    genetic_matrix = np.zeros((n_snp,n_pop)) #allele freq for each population
    
    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    #Population matrix
    pop_matrix = np.zeros((n_pop,n_indiv)) #individual population admixture
    
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
        genetic_matrix[i,:] = np.random.beta(frq[i]*(1-fst[i])/fst[i], (1-frq[i])*(1-fst[i])/fst[i], size=n_pop)            
        
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
    return (pop_matrix, pop_idx, genetic_matrix)

def model_psd(n_pop, n_snp, n_indiv, frq, fst, seed: None):
    '''Simulate model psd.
        Parameters
        ----------
        n_pop : int
            Number of population(s).
        n_snp : int
            Number of SNPs. 
        n_indiv : int
            Number of individuals.
        frq : array
            Frequency/proportion of a specific SNP in the data
        fst : array
            Percentage of contribution to the total generic variation from each SNP
        seed : int
            Uses datetime if nothing is provided
            Randomizes generation.
        Returns
        -------
        pop_matrix : 2d data array
            Individual population adimixture matrix.
        pop_idx : array
            Population adimixture matrix.
        genetic_matrix : 2d data array
            Genetic data matrix.
    '''

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %%%%%%% Load HapMap Info
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #HM_inf = pd.read_csv(u'CEUASWMEX_fst_frq.txt',sep=u' ')
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
    #d = 3 #number of populations
    # % each row of Gamma will be populated with 3 i.i.d. draws from BN model
    # % thanks to the nature of the HapMap data (sampled from 3 populations)
    #
    # % allele freq of population: allele freq of each SNP described by that
    # % population
    genetic_matrix = np.zeros((n_snp,n_pop)) #allele freq for each population

    # % columns of S will be populated with indicator vectors s.t. each
    # % individual assigned to one of the 3 subpopulations i.e. admixture
    pop_matrix = np.zeros((n_pop,n_indiv)) #individual population admixture

    # %X = zeros(m,n); %genotype matrix
    if (seed == None):
        seed = datetime.now()
    
    # % random seeding here...
    random.seed(seed)

    # %populate the allele freq matrix from BN with (p,F) from HapMap
    # % for each SNP...
    for i in range(0,n_snp):
        # each row will generate 'd' variates drawing from this distribution
        genetic_matrix[i,:] = np.random.beta(frq[i]*(1-fst[i])/fst[i], (1-frq[i])*(1-fst[i])/fst[i], size=n_pop)

    # print('Mean of allele freq matrix: ', np.mean(G, axis=0))

    # %populate the population admixture matrix
    # %this part is tailor made for 3 populations as described in
    # %Song et al. Nat. Genet. (2015)
    # % Treating the probabilities as ranges
    # % 1: <60/210, 2: bet. 60/210 and 120/210, 3: >=120/210

    alpha = 0.1*np.ones((n_pop,1))
    popidx = np.zeros((n_indiv,1))
    for i in range(0,n_indiv):
        for j in range(0,n_pop):
            pop_matrix[j,i] = np.random.gamma(alpha[j],1)

        pop_matrix[:,i] = pop_matrix[:,i]/np.sum(pop_matrix[:,i])
        I = np.argmax(pop_matrix[:,i])
        popidx[i] = I+1

    return (pop_matrix, popidx, genetic_matrix)

def model_HGDP(flag, n_pop, n_snp, n_indiv, topPCs, HGDP_subpops, seed: None):
    '''Simulate model HGDP
        Parameters
        ----------
        flag : int
            Plot flag.
        n_pop : int
            Number of populations.
        n_snp : int
            Number of SNPs
        n_indiv : int
            Number of individuals.
        topPCs : 2d data array
            TBD
        HGDP_subpops : array
            TBD
        seed : int
            Uses datetime if nothing is provided
            Randomizes generation.
        
        Returns
        -------
        pop_matrix : 2d data array
            Individual population adimixture matrix.
        pop_idx : array
            Population adimixture matrix.
        genetic_matrix : 2d data array
            Genetic data matrix.
    '''

    # REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # Downsampling for simulation (computationally easier)
    # m = int(1e4) #number of SNPs
    # n = int(305) #no. of individuals
    #flag = 0 #plot flag
    #d = int(10) #number of populations (see log of --fst result)

    # allele freq of population: allele freq of each SNP described by that
    # population
    genetic_matrix = np.zeros((n_snp,n_pop)) #allele freq for each population

    # columns of S will be populated with indicator vectors s.t. each
    # individual assigned to one of the 51 subpopulations i.e. admixture
    pop_matrix = np.zeros((n_pop,n_indiv)) #individual population admixture

    # Uses current time if no seed is provided
    if (seed == None):
        seed = datetime.now()
    
    # % random seeding here...
    random.seed(seed)

    #populate the allele freq matrix from BN with (p,F) from HapMap
    # for each SNP...
    for i in range(0,n_snp):
    # each row will generate 'd' variates drawing from this distribution
        genetic_matrix[i,:] = 0.9*np.random.uniform(0, 0.5, size=n_pop)

    # set last column to 0.05 per Song et al. 2015
    genetic_matrix[:,n_pop-1] = 0.05

    for i in range(0,n_pop):
        pop_matrix[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
    pop_matrix[n_pop-1,:] = 1

    popidx = np.zeros((n_indiv,1))

    # HGDP_subpops = pd.read_csv(u'subpops_pruned_HGDP.txt',sep=u' ',header=None)

    for i in range(0,HGDP_subpops.shape[0]):
        if HGDP_subpops[i] == u"Biaka_Pygmies":
            popidx[i] = 1
        elif HGDP_subpops[i] == u"French":
            popidx[i] = 2
        elif HGDP_subpops[i] == u"Han":
            popidx[i] = 3
        elif HGDP_subpops[i] == u"Japanese":
            popidx[i] = 4
        elif HGDP_subpops[i] == u"Palestinian":
            popidx[i] = 5
        elif HGDP_subpops[i] == u"Papuan":
            popidx[i] = 6
        elif HGDP_subpops[i] == u"Pima":
            popidx[i] = 7
        elif HGDP_subpops[i] == u"Russian":
            popidx[i] = 8
        elif HGDP_subpops[i] == u"Sardinian":
            popidx[i] = 9
        else:
            # Sindhi
            popidx[i] = 10

    return (pop_matrix, popidx, genetic_matrix)

def model_TGP(flag, n_pop, n_snp, n_indiv, topPCs, TGP_subpops, __, seed: None):
    '''Simulate model tgp.
        Parameters
        ----------
        flag : int
            Plot flag.
        n_pop : int
            Number of populations.
        n_snp : int
            Number of SNPs
        n_indiv : int
            Number of individuals.
        topPCs : 2d data array
            TBD
        TGP_subpops : array
            TBD
        seed : int
            Uses datetime if nothing is provided
            Randomizes generation.
        Returns
        -------
        pop_matrix : 2d data array
            Individual population adimixture matrix.
        pop_idx : array
            Population adimixture matrix.
        genetic_matrix : 2d data array
            Genetic data matrix.
    '''
    # REMEMBER: 'n' here is INDIVIDUALS not SNPs
    # Downsampling for simulation (computationally easier)
    # m = int(1e4) #number of SNPs
    # n = int(1056) #no. of individuals
    # flag = 0 #plot flag
    # d = int(10) #number of populations (see log of --fst result)

    # allele freq of population: allele freq of each SNP described by that
    # population
    genetic_matrix = np.zeros((n_snp,n_pop)) #allele freq for each population

    # columns of S will be populated with indicator vectors s.t. each
    # individual assigned to one of the 51 subpopulations i.e. admixture
    pop_matrix = np.zeros((n_pop,n_indiv)) #individual population admixture

    # Uses current time if no seed is provided
    if (seed == None):
        seed = datetime.now()
    
    # % random seeding here...
    random.seed(seed)

    #populate the allele freq matrix from BN with (p,F) from HapMap
    # for each SNP...
    for i in range(0,n_snp):
        # each row will generate 'd' variates drawing from this distribution
        genetic_matrix[i,:] = 0.9*np.random.uniform(0, 0.5, size=n_pop)

    # set last column to 0.05 per Song et al. 2015
    genetic_matrix[:,n_pop-1] = 0.05

    # TGP_PCs = pd.read_csv(u'pruned_TGP_topPops_singVecs.txt',sep=u' ',header=None)
    # topPCs = TGP_PCs.values

    for i in range(0,n_pop):
        pop_matrix[i,:] = (topPCs[:,i]-np.min(topPCs[:,i]))/(np.max(topPCs[:,i])-np.min(topPCs[:,i]))
    pop_matrix[n_pop-1,:] = 1

    popidx = np.zeros((n_indiv,1))

    # TGP_subpops = pd.read_csv(u'subpops_pruned_TGP.txt',sep=u' ',header=None)

    for i in range(0,TGP_subpops.shape[0]):
        if TGP_subpops[i] == u"CHB":
            popidx[i] = 1
        elif TGP_subpops[i] == u"CHS":
            popidx[i] = 2
        elif TGP_subpops[i] == u"GIH":
            popidx[i] = 3
        elif TGP_subpops[i] == u"GWD":
            popidx[i] = 4
        elif TGP_subpops[i] == u"IBS":
            popidx[i] = 5
        elif TGP_subpops[i] == u"JPT":
            popidx[i] = 6
        elif TGP_subpops[i] == u"PUR":
            popidx[i] = 7
        elif TGP_subpops[i] == u"STU":
            popidx[i] = 8
        elif TGP_subpops[i] == u"TSI":
            popidx[i] = 9
        else:
            # YRI
            popidx[i] = 10

    return (pop_matrix, popidx, genetic_matrix)
