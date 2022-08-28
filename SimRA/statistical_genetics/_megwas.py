import numpy as np


def generate_genotypes(model_flag, genetic_matrix, pop_matrix):
    """generate_genotypes 
        Parameters
        ----------
        model_flag : str
            Model flag of population.
        genetic_matrix : 2d data array
            Genetic data matrix.
        pop_matrix : 2d data array
            Individual population admixture matrix.
        
        Returns
        -------
        n_pop : int
            Number of populations.
        geno_matrix : 2d data array
            Output genotype matrix.
            
    """
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
    return n_pop, geno_matrix