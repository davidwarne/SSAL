
#include "SSAL_sequential.h"
#include "SSAL_util.h"
#include "SSAL_RNG.h"
/**
 * @brief Single Precision Gillespie exact stochastic simulation algorithm (eSSA)
 * @details Simulates a discrete-state continuous time Markov process.
 *
 * @param m number of reactions
 * @param n dimension of the state vector
 * @param nt the number of timesteps to sample at
 * @param T an array of length nt sample times
 * @param X0 initial state vector 
 * @param nu_minus reactant stoichiometric coefficients m x n (stored in row-major format)
 * @param nu the stoichiometric matrix m x n ( store in row-major format)
 * @param c kinetic reaction rates
 * @param ndims dimension to measure
 * @param dims array of indices of length ndims
 * @param X_r state-space trajection for measured dims (nt*ndims)
 *
 * @retVal 
 */
int segils(int m,int n,int nt,float * restrict T, float * restrict X0, float *restrict nu_minus,
    float * restrict nu,float * restrict c,int ndims,int *restrict dims,float *restrict X_r)
{
    float a[m]; /*propensities or hazard values*/
    float X[n]; /* state vector*/
    float t=0;
    float a_0; /*combined propensity*/
    float deltat; /* time to next reaction*/
    float r1,r2; /*random numbers*/

    int i,j,ti;

    /*copy initial condition*/
    for (i=0;i<n;i++)
    {
        X[i] = X0[i];
    }

    /*compute propensities*/
    suhzds(m,n,nu_minus,X,a);

    /*total propensity*/
    a_0 = 0;
    for (j=0;j<n;j++)
    {
        a_0 += a[j];
    }

    /*generate time to next reaction dt ~ Exp(a_0)*/
    r1 = ((float)rand())/((float)ONE_ON_RAND_MAX);
    deltat = -log(r1)/a_0;

    t=0;
    for (ti=0;ti<nt;ti++)
    {
        /*run a gillespie till we reach the next measure time*/
        while (t <= T[ti])
        {
            float a_sum;
            int k;
            /* sample discrete distribution to determine the reaction
             * channel 
             */
            r2 = ((float)rand())/((float)ONE_ON_RAND_MAX);
            a_sum = 0;
            k = 0;
            while (a_sum <= a_0*r2 && k < m - 1)
            {
                k++
                a_sum += a[k]
            }

            /*update state*/
            for (i=0;i<n;i++)
            {
                X[i] += nu[n*k + i];
            }
            /*update time*/
            t = t + deltat;

            /*update propensities*/
            suhzds(m,n,nu_minus,c,X,a);
            /*total propensity*/
            a_0 = 0;
            for (j=0;j<n;j++)
            {
                a_0 += a[j];
            }
        
            /*generate time to next reaction dt ~ Exp(a_0)*/
            r1 = ((float)rand())/((float)ONE_ON_RAND_MAX);
            deltat = -log(r1)/a_0;
        }

        /*record our next measurement*/
        for (i=0;i<ndims;i++)
        {
            X_r[ndims*ti + i] = X[dims[i]];
        }
    } 

    return 0;
}
