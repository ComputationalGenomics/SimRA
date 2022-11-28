import numpy as np

import _genotypes


def generate_genotypes(model_flag, n_snp, n_indiv, topPCs: None, subpops: None, flag: None, seed: None):
    """generate_genotypes 
        Parameters
        ----------
        model_flag : str
            Model flag of population.
        n_snp : int
            Number of SNPs
        n_indiv : int
            Number of individuals.
        topPCs : 2d data array
            Only for HGDP & TGP - TBD.
        subpops : array
            Only for HGDP & TGP - TBD.
        flag : int
            Only for HGDP & TGP - Plot flag.
        seed : int
            Optional - Randomizes generation
            Uses datetime if nothing is provided.
        Returns
        -------
        geno_matrix : 2d data array
            Output genotype matrix.
            
    """
    if model_flag == u'BN':
        pop_matrix, popidx, genetic_matrix = _genotypes.model_bn(3, )
    elif model_flag == u'PSD':
        pop_matrix, popidx, genetic_matrix = _genotypes.model_psd(3, )
    elif model_flag == u'HGDP':
        pop_matrix, popidx, genetic_matrix = _genotypes.model_hgdp(flag, 10)
    elif model_flag == u'TGP':
        pop_matrix, popidx, genetic_matrix = _genotypes.model_tgp(flag, 10)

    # %Get the allele frequency matrix for each individual (SNP -by- indiv)
    # % GOAL: efficiently estimate the individual-specific frequencies, F
    # % Hao 2015 paper
    allele_freq = np.matmul(genetic_matrix,pop_matrix)

    if model_flag == u"HGDP" or u"TGP":
        #####################################
        # Normalize allele_freq (formerly F) by column (making sure each column i.e. individual is bet. [0 1])
        allele_freq = allele_freq/allele_freq.max(axis=0)
        n_pop = 10
    else:
        n_pop = 3
        #####################################

    # % simulating X using binomial distribution of estimated F
    geno_matrix = np.random.binomial(2, allele_freq)

    # # % if A is a matrix, then sum(A,2) is a column vector containing the sum of each row.
    idxzer = np.where(~geno_matrix.any(axis=1))[0]
    # print(idxzer)
    return geno_matrix