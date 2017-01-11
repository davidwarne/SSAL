/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2016  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ASSA_sequential.h"
#include "util_sequential.h"
#include <time.h>

/**
 * @brief Double Precision Multi-level configures level sample sizes
 * @details uses a predefined fixed number of samples for each level l to
 * determine the optimal number of samples per level to minimise runtime 
 * given a target estimator variance.
 *
 * @param m number of reactions
 * @param n state space vector dimension
 * @param T simulation endtime T
 * @param X0 initial conditions
 * @param nu_minus reactant stoichiometric coefficients
 * @param nu the stoichiometric matrix 
 * @param c kinetic reaction rates
 * @param tau0 base level tau (coare level)
 * @param M scale factor tau_l = tau0*M^(-l)
 * @param L number of levels
 * @param epsilon target estimator standard deviation
 * @param ndims the number of state dimensions to measure
 * @param dim array of state dimensions indexes
 * @param f function of state vector f : \mathbb{Z}^n -> \mathbb{R}^n
 * @param nl array of length L of optimal sample numbers
 */
int dumlnls(int m, int n, double T, double * restrict X0, double * restrict nu_minus,
    double * restrict nu, double * restrict c, double tau0, int M, int L, double epsilon,
    int ndims, int * restrict dims, int (*f)(int, double *, double *), int * restrict nl)
{
    double n_f = 0;
    double E_l[ndims];
    double E_l2[ndims];
    double taul;
    double Z_l_r[ndims];
    double Z_lm1_r[ndims];
    double fZ_l_r[ndims];
    double fZ_lm1_r[ndims];
    double times[L];
    double sigma2[L];
    int i,j,l;
    taul = tau0;
    /*time N trials of level-l and compute variances*/
    for (l=0;l<L;l++)
    {
        clock_t start_t;

        for (i=0;i<ndims;i++)
        {
            E_l[i] = 0;
        }
        for (i=0;i<ndims;i++)
        {
            E_l2[i] = 0;
        }

        start_t = clock();        
        if (f == NULL) /* use E[Z(T)] instead of E[f(Z(T))]*/
        {
            if (l != 0)
            {
                for (j=0;j<MLMC_TIMING_TRIALS;j++)
                {
                    dactauls(m,n,1,&T,X0,nu_minus,nu,c,ndims,dims,taul,M,Z_l_r,Z_lm1_r);
                    for (i=0;i<ndims;i++)
                    {
                        E_l[i] += (Z_l_r[i] - Z_lm1_r[i]);
                    }

                    for (i=0;i<ndims;i++)
                    {
                        E_l2[i] += (Z_l_r[i] - Z_lm1_r[i])*(Z_l_r[i] - Z_lm1_r[i]);
                    }
                }
            }
            else
            {
                for (j=0;j<MLMC_TIMING_TRIALS;j++)
                {
                    datauls(m,n,1,&T,X0,nu_minus,nu,c,ndims,dims,taul,Z_l_r);
                    for (i=0;i<ndims;i++)
                    {
                        E_l[i] += Z_l_r[i];
                    }

                    for (i=0;i<ndims;i++)
                    {
                        E_l2[i] += Z_l_r[i]*Z_l_r[i];
                    }
                }
            }

        }
        else
        {
            if (l != 0)
            {
                for (j=0;j<MLMC_TIMING_TRIALS;j++)
                {
                    dactauls(m,n,1,&T,X0,nu_minus,nu,c,ndims,dims,taul,M,Z_l_r,Z_lm1_r);
                    (*f)(ndims,Z_l_r,fZ_l_r);
                    (*f)(ndims,Z_lm1_r,fZ_lm1_r);
                    for (i=0;i<ndims;i++)
                    {
                        E_l[i] += (fZ_l_r[i] - fZ_lm1_r[i]);
                    }

                    for (i=0;i<ndims;i++)
                    {
                        E_l2[i] += (fZ_l_r[i] - fZ_lm1_r[i])*(fZ_l_r[i] - fZ_lm1_r[i]);
                    }
                }

            }
            else
            {
                for (j=0;j<MLMC_TIMING_TRIALS;j++)
                {
                    datauls(m,n,1,&T,X0,nu_minus,nu,c,ndims,dims,taul,Z_l_r);
                    (*f)(ndims,Z_l_r,fZ_l_r);
                    for (i=0;i<ndims;i++)
                    {
                        E_l[i] += fZ_l_r[i];
                    }

                    for (i=0;i<ndims;i++)
                    {
                        E_l2[i] += fZ_l_r[i]*fZ_l_r[i];
                    }
                }
            }
        }
        /*not the best of timers, but should be sufficient*/
        times[l] = (double)(clock() - start_t)/((double)CLOCKS_PER_SEC); 
        times[l] /= (double)MLMC_TIMING_TRIALS;        
        /*compute marginal sample sigma^2*/
        for (i=0;i<ndims;i++)
        {
            E_l2[i] /= (double)MLMC_TIMING_TRIALS;
        }
        for (i=0;i<ndims;i++)
        {
            E_l[i] /= (double)MLMC_TIMING_TRIALS;
        }
        for (i=0;i<ndims;i++)
        {
            E_l2[i] -= E_l[i]*E_l[i];
        }
        /*the sigma^2 we use is the dim with max sigma^2*/
        sigma2[l] = E_l2[0];

        for (i=1;i<ndims;i++)
        {
            if (E_l2[i] > sigma2[l])
            {
                sigma2[l] =  E_l2[i];
            }
        }
        
        taul /= (double)M;
    }

    /*now we can compute "optimal nl" derived from Legrange multiplier
     * e.g., Minimize f(n_0,...,n_l) = \sum_{l=0}^L c_l n_l
     *       Subject to g(n_0,...,n_l) = \sum_{l=0}^L \sigma_l^2/n_l = \epsilon^2
     * where c_l is the compute time of a realisation of Z_l(T) and \sigma_l^2 is
     * the sample variance of Z_l(T).
     *
     * the solution is,
     *      n_l = \sigma_l / (\epsilon^2 \sqrt(c_l)) \sum_{m=0}^L \sigma_m \sqrt(c_m)
     */
    for (l=0;l<L;l++)
    {
        n_f += sqrt(sigma2[l]*times[l]);
    }

    for (l=0;l<L;l++)
    {
        nl[l] = (int)ceil(n_f*(sqrtf(sigma2[l])/sqrtf(times[l]))/(epsilon*epsilon));
    }

    return 0;
}
