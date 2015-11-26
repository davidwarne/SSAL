/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2015  David J. Warne
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
 * @brief Single precision Multi-level Monte Carlo method
 * @datails An efficient Monte Carlo estimator with the same statistical 
 * bias against exact simulation as a tau-leaping with fine-grain tau.
 * This function computes for a tau leap approximation Z_L with tau = tau0*M^(-L),
 *      E[f(Z_L(T))] ~ Q_0 + sum_l=1^L Q_l + O(epsilon)
 * where
 *      Q_0 = 1/n0 sum_j=1^n0 f(Z_0^j(T)) with Z_0^j the j-th realisation of Z_0
 *      Q_l = 1/nl sum_j=1^nl f(Z_l^j(T)) - f(Z_{l-1}^j(T))
 *
 * The function is an unbiased esitmator for E[f(Z_L(T))] thus,
 *      E[f(X(T))] = E[f(Z_L(T))] + O(tauL)
 *
 * Using the desired estimator error epsilon a small number of samples are utilised 
 * determine optimal choices for the sample numbers
 *
 * @param m number of reactions
 * @param n dimension of state vector
 * @param nt number of timesteps to sample at
 * @param T arr of length nt of sample times
 * @param X0 Initial condition
 * @param nu_minus reactant stoichiometric coefficients m x n (row-major format)
 * @param nu the stoichiometric matrix m x n (row-major format)
 * @param c the kinetic reaction rates
 * @param tau0 base level tau such that taul = tau0*M^(-l)
 * @param M scale factor between levels
 * @param L the number of levels
 * @param epsilon the desired statistical error (desired standard deviation)
 * @param ndims numebr dimension to measure
 * @param dims the dimensions indices 
 * @param f functional of state vector f : \mathbb{Z}^n -> \mathbb{R}^n 
 * @param E_X expected value of f(X(T))
 * @param V_X varience of f(X(T))
 *
 * @param if f == NULL then the expected state vector is used
 */
int samlmcbs(int m,int n, int nt, float * restrict T, float * restrict X0, float * restrict nu_minus,
    float * restrict nu, float * restrict c, float tau0, int M, int L, float epsilon, 
    int ndims,int * restrict dims, int (*f)(int,float *, float *), float * restrict E_X, float * restrict V_X)
{
    int i,j,l; 
    int nl[L]; /*sample sizes for each level*/
    float taul;
    float times[L];
    float sigma2[L];
    float E_l[nt*ndims];
    float E_l2[nt*ndims];
    float Z_l_r[nt*ndims];
    float fZ_l_r[nt*ndims];
    float Z_lm1_r[nt*ndims];
    float fZ_lm1_r[nt*ndims];

    for (i=0;i<nt*ndims;i++)
    {
        E_X[i]  = 0;
    }
    for (i=0;i<nt*ndims;i++)
    {
        V_X[i] = 0;
    }

    /*timing analysis for each level*/
    sumlnls(m,n,T[nt-1],X0,nu_minus,nu,c,tau0,M,L,epsilon,ndims,dims,f,nl);
    if (f == NULL) /* use E[Z(T)] instead of E[f(Z(T))]*/
    {
        /*run level 0 samples*/
        for (j=0;j<nl[0];j++)
        {
            /*Tau-leaping for Z_0(t) for t = T[0],...,T[nt-1]*/
            satauls(m,n,nt,T,X0,nu_minus,nu,c,ndims,dims,tau0,Z_l_r);
            /*accumulate Z*/
            for (i=0;i<nt*ndims;i++)
            {
                E_X[i] += Z_l_r[i];
            }

            for (i=0;i<nt*ndims;i++)
            {
                V_X[i] += Z_l_r[i]*Z_l_r[i];
            }
        }
    
        /*level 0 approximation with high statistical biased*/
        for (i=0;i<nt*ndims;i++)
        {
            E_X[i] /= (float)nl[0];
        }
        /*variance of level 0*/
        for (i=0;i<nt*ndims;i++)
        {
            V_X[i] /= (float)nl[0];
        }
        for (i=0;i<nt*ndims;i++)
        {
            V_X[i] -= E_X[i]*E_X[i];
        }
        for (i=0;i<nt*ndims;i++)
        {
            V_X[i] /= (float)nl[0];
        }
        
    
        taul = tau0/((float)M);
        /*run each bias correction term*/
        for (l=1;l<L;l++)
        {
            if (nl[l] == 0)
            {
                taul /= (float)M;
                continue;
            }
            for (i=0;i<nt*ndims;i++)
            {
                E_l[i] = 0;
            }
            for (i=0;i<nt*ndims;i++)
            {
                E_l2[i] = 0;
            }
            for (j=0;j<nl[l];j++)
            {
                /*run correlated  Z_l(t),Z_{l-1}(t)*/ 
                sactauls(m,n,nt,T,X0,nu_minus,nu,c,ndims,dims,taul,M,Z_l_r,Z_lm1_r);
                for (i=0;i<nt*ndims;i++)
                {
                    E_l[i] += (Z_l_r[i] - Z_lm1_r[i]);
                }
                for (i=0;i<nt*ndims;i++)
                {
                    E_l2[i] += (Z_l_r[i] - Z_lm1_r[i])*(Z_l_r[i] - Z_lm1_r[i]);
                }
            }
    
            /*level l bias correction*/
            for (i=0;i<nt*ndims;i++)
            {
                E_l[i] /= (float)nl[l];
            }
            /*level l corrector variance*/
            for (i=0;i<nt*ndims;i++)
            {
                E_l2[i] /= (float)nl[l];
            }
            for (i=0;i<nt*ndims;i++)
            {
                E_l2[i] -= E_l[i]*E_l[i];
            }
            /*add to E_X*/
            for (i=0;i<nt*ndims;i++)
            {
                E_X[i] += E_l[i];
            }

            /*add to V_X*/
            for (i=0;i<nt*ndims;i++)
            {
//                V_X[i] += E_l2[i]/((float)nl[l]);
            }
            taul /= (float)M;
        }

    }
    else
    {
        /*run level 0 samples*/
        for (j=0;j<nl[0];j++)
        {
            /*Tau-leaping for Z_0(t) for t = T[0],...,T[nt-1]*/
            satauls(m,n,nt,T,X0,nu_minus,nu,c,ndims,dims,tau0,Z_l_r);
            /*evaluate functional f(Z_0(t))*/
            (*f)(n,Z_l_r,fZ_l_r);
            /*accumulate fZ*/
            for (i=0;i<nt*ndims;i++)
            {
                E_X[i] += fZ_l_r[i];
            }
        }
    
        /*level 0 approximation with high statistical biased*/
        for (i=0;i<nt*ndims;i++)
        {
            E_X[i] /= (float)nl[0];
        }
    
        taul = tau0/((float)M);
        /*run each bias correction term*/
        for (l=1;l<L;l++)
        {
            for (i=0;i<nt*ndims;i++)
            {
                E_l[i] = 0;
            }
            for (j=0;j<nl[l];j++)
            {
                /*run correlated  Z_l(t),Z_{l-1}(t)*/ 
                sactauls(m,n,nt,T,X0,nu_minus,nu,c,ndims,dims,taul,M,Z_l_r,Z_lm1_r);
                /*evaluate functional f(Z_l(t) and f(Z_{l-1}(t))*/
                (*f)(n,Z_l_r,fZ_l_r);
                (*f)(n,Z_lm1_r,fZ_lm1_r);
                for (i=0;i<nt*ndims;i++)
                {
                    E_l[i] += (fZ_l_r[i] - fZ_lm1_r[i]);
                }
                for (i=0;i<nt*ndims;i++)
                {
                    E_l2[i] += (fZ_l_r[i] - fZ_lm1_r[i])*(fZ_l_r[i] - fZ_lm1_r[i]);
                }
            }
    
            /*level l bias correction*/
            for (i=0;i<nt*ndims;i++)
            {
                E_l[i] /= (float)nl[l];
            }
            /*level l corrector variance*/
            for (i=0;i<nt*ndims;i++)
            {
                E_l2[i] /= (float)nl[l];
            }
            for (i=0;i<nt*ndims;i++)
            {
                E_l2[i] -= E_l[i]*E_l[i];
            }
            /*add to E_X*/
            for (i=0;i<nt*ndims;i++)
            {
                E_X[i] += E_l[i];
            }
            /*add to V_X*/
            for (i=0;i<nt*ndims;i++)
            {
                V_X[i] += E_l2[i]/((float)nl[l]);
            }
            taul /= (float)M;
        }
    }
    return 0;
}
