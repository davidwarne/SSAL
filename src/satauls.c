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

/**
 * @detail Single precision tau-leap method
 * @detail An approximate stochastic simulation algorithm in which a continuous time
 * markov process is discretised in time
 *
 * @param m number of reactions
 * @param n dimension of state vector
 * @param nt the number of timesteps to sample at
 * @param T an array of length nt of sample times
 * @param X0 Initial condition.
 * @param nu_minus reactant stoichiometric coefficients m x n (row-major format)
 * @param nu the stoichiometric matrix m x n (row-major format)
 * @param c the kinetic reaction rates
 * @param ndims numebr dimension to measure
 * @param dims the dimensions indices 
 * @param X_r state-space trajectory for measured dims (nt*ndims)
 * @param tau the discrete timestep parameter
 */
int satauls(int m,int n,int nt,float * restrict T, float * restrict X0, float *restrict nu_minus,
    float * restrict nu,float * restrict c,int ndims,int *restrict dims,float *restrict X_r,
    float tau)
{
    float a[m]; /*propensities*/
    float X[n]; /* the state vector*/
    float t=0;
    float p[m]; /* poisson random variates */
    int i,j,ti;

    for (i=0;i<n;i++)
    {
        X[i] = X0[i];
    }

    suhzds(m,n,nu_minus,c,X,a);

    for (ti=0;ti<nt;ti++)
    {
        for (;t <= T[ti];t+=tau)
        {
            /*generate poisson variates*/
            for (j=0;j<m;j++)
            {
                p[j] = (float)surngpois(a[j]*tau);
         
            }
            /*update state vector*/
            for (j=0;j<m;j++)
            {
                for (i=0;i<n;i++)
                {
                    X[i] += p[j]*nu[j*n +i];
                }
            }
            /*update propensities*/
            suhzds(m,n,nu_minus,c,X,a);
        }

        /*write out timestep*/
        for (i=0;i<ndims;i++)
        {
            X_r[i*nt+ti] = X[dims[i]];
        }
    }

}

