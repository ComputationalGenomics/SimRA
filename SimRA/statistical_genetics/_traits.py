import numpy as np
import math
from sklearn.cluster import KMeans

################################################################################
################################################################################
################################################################################

def gen_eff (gtype_matrix, variance_gen, num_indiv, num_causal_snp=10) :
    """Simulates genetic effect on each individual's quantitative traits and disease susceptibility.

        Parameters
        ----------
        gtype_matrix : 2-d data array
            Normalized genotype matrix.
        variance_gen : int
            Amount of variance for genetic effect. 
            ***When simulating for all 3 effects, variances for each must sum to 100***
        num_indiv : int
            Number of individuals.
        num_causal_snp : int, default 10
            Number of causal SNPs. Default is 10.

        Returns
        -------
        gen_eff : 1-d data array
            Simulated genetic effect on each individual's quantitative traits and disease susceptibility.

        See Also
        --------
        enviro_eff
        noise_eff
    """

    variance_gen /= 100.0

    snp_effects = np.random.normal(0,0.5,num_causal_snp)

    # genetic effect
    for i in range (0,num_causal_snp) :
        gtype_matrix[i,:] *= snp_effects[i]

    # sum of each column of normX
    gen_eff = gtype_matrix.sum(axis = 0)

    # rescale genetic effects by the variance factor
    fact = math.sqrt(variance_gen)/(np.std(gen_eff))
    # element-wise rescaling of the genetic effects
    gen_eff *= fact

    gen_eff = gen_eff.reshape((num_indiv,1))

    return gen_eff
    


def enviro_eff (indiv_pop_admixture, variance_enviro, num_indiv, num_pop) :
    """Simulates environmental effect on each individual's quantitative traits and disease susceptibility.

        Parameters
        ----------
        indiv_pop_admixture : 2-d data array
            Individual population admixture matrix.
        variance_enviro : int
          Amount of variance for environmental effect.
            ***When simulating for all 3 effects, variances for each must sum to 100***
        num_indiv : int
            Number of individuals.
        num_pop : int
            Number of population(s).

        Returns
        -------
        enviro_eff : 2-d data array
            Simulated environmental effect on each individual's quantitative traits and disease susceptibility.

        See Also
        --------
        gen_eff
        noise_eff
    """

    variance_enviro /= 100.0

    # k-means clustering of the columns of the population admixture matrix
    # indiv_pop_admixture.T each individual falls in one of 'num_pop' clusters
    theclustering = KMeans(n_clusters=num_pop).fit(indiv_pop_admixture.T)
    # print(theclustering.labels_)
    enviro_eff = theclustering.labels_ + 1
    

    # rescale enviro_eff
    mc_enviro_eff = np.square((enviro_eff - np.mean(enviro_eff)))
    fact = math.sqrt(variance_enviro)/math.sqrt(np.sum(mc_enviro_eff)/(num_indiv-1))
    # element-wise rescaling of the environmental effects
    enviro_eff *= fact

    enviro_eff = enviro_eff.reshape((num_indiv,1))

    return enviro_eff


def noise_eff (indiv_pop_admixture, variance_noise, num_indiv, num_pop) :
    """Simulates noise effect on each individual's quantitative traits and disease susceptibility.

        Parameters
        ----------
        indiv_pop_admixture : 2-d data array
            Individual population admixture matrix.
        variance_noise : int
          Amount of variance for noise effect.
            ***When simulating for all 3 effects, variances for each must sum to 100***
        num_indiv : int
            Number of individuals.
        num_pop : int
            Number of population(s).

        Returns
        -------
        noise_eff : 2-d data array
            Simulated noise effect on each individual's quantitative traits and disease susceptibility.

        See Also
        --------
        gen_eff
        enviro_eff
    """

    variance_noise /= 100.0

    noise_eff = np.zeros((num_indiv,1))

    # k-means clustering of the columns of the population admixture matrix
    # s.t. each individual falls in one of 'd' clusters
    # combine the 2 lines or make individual? which is better for readability?
    theclustering = KMeans(n_clusters=num_pop).fit(indiv_pop_admixture.T)
    # print(theclustering.labels_)

    # some gamma distribution for the noise
    gamvec = 1/np.random.gamma(3,1,size=(num_pop,1))

    for i in range(0,num_indiv):
        noise_eff[i] = np.random.normal(0,gamvec[theclustering.labels_[i]-1])

    # partly replaceable with np.std?
    mc_noise_eff = np.square((noise_eff - np.mean(noise_eff)))
    fact = math.sqrt(variance_noise)/math.sqrt(np.sum(mc_noise_eff)/(num_indiv-1))
    # element-wise rescaling of the noise effects
    noise_eff *= fact

    return noise_eff


def combined_effs (gtype_matrix, indiv_pop_admixture, all_variances, num_indiv, num_pop, num_causal_snp=10) :
    """Simulates quantitative traits and disease susceptibility for each individual.

        Parameters
        ----------
        gtype_matrix : 2-d data array
            Normalized genotype matrix.
        indiv_pop_admixture : 2-d data array
            Individual population admixture matrix.
        all_variances : int tuple
          Amount of variance for each effect (genetic, environmental and noise in that order).
            ***Variances for each effect must sum to 100***
        num_indiv : int
            Number of individuals.
        num_pop : int
            Number of population(s).
        num_causal_snp : int, default 10
            Number of causal SNPs. Default is 10.

        Returns
        -------
        traits : vector
            Simulated quantitative trait for each individual.
        status : vector
            Simulated disease susceptibility (binary trait) for each individual (case status (1) or control (0))

        See Also
        --------
        gen_eff
        enviro_eff
        noise_eff
    """

    # should probably raise an exception/do try: except: for the variance not summing to 100 (also is it a float or int?)

    gen_eff = gen_eff (gtype_matrix, all_variances[0], num_causal_snp, num_indiv)

    enviro_eff = enviro_eff (indiv_pop_admixture, all_variances[1], num_indiv, num_pop)

    noise_eff = noise_eff (indiv_pop_admixture, all_variances[2], num_indiv, num_pop)

    # Get the simulated quantitative trait for each individual
    # apostrophe --> transpose
    traits = gen_eff + enviro_eff + noise_eff
    print(' ')

    # Binary traits
    bin_trait = traits-noise_eff
    # probability that the individual has case status (vs control)
    prob_case_status = np.power(10,bin_trait)/(1+np.power(10,bin_trait))
    # print(np.where(prob_case>0.5)[0])
    status = np.zeros((num_indiv,1))
    status[np.where(prob_case_status > 0.5)[0]] = 1

    return traits, status


if __name__ == '__main__':
    print("HEY HEY HEY")