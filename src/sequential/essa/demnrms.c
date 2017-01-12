/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2017  David J. Warne
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
#include "ESSA_sequential.h"
#include "util_sequential.h"

#include <stdio.h>

/**
 * @brief Double Precision Modified Next Reaction Method  exact stochastic 
 * simulation algorithm(eSSA)
 * @details Simulates a discrete-state continuous time Markov process.
 *
 * @param m number of reactions
 * @param n dimension of the state vector
 * @param nt the number of timesteps to sample at
 * @param T an array of length nt sample times
 * @param X0 initial state vector 
 * @param nu_minus reactant stoichiometric coefficients m x n 
 *                 (stored in row-major format)
 * @param nu the stoichiometric matrix m x n ( store in row-major format)
 * @param c kinetic reaction rates
 * @param ndims dimension to measure
 * @param dims array of indices of length ndims
 * @param X_r state-space trajection for measured dims (nt*ndims)
 *
 * @retVal 
 */
int 
demnrms(int m,int n,int nt,double * restrict T, double * restrict X0, 
       double *restrict nu_minus, double * restrict nu,double * restrict c,
       int ndims,int *restrict dims,double *restrict X_r)
{
    double a[m]; /*propensities or hazard values*/
    double X[n]; /* state vector*/
    double t=0;
    double delta;
    double deltat[m]; /* time to next reaction*/
    double P[m];
    double T_r[m];
    int mu;
    int i,j,ti;

    /*copy initial condition*/
    for (i=0;i<n;i++)
    {
        X[i] = X0[i];
    }

    for (j=0;j<m;j++)
    {
        T_r[j] = 0;
    }

    /*compute propensities*/
    duhzds(m,n,nu_minus,c,X,a);
    

    for (j=0;j<m;j++)
    {
        P[j] = durngexps(1);
    }
    for (j=0;j<m;j++)
    {
        deltat[j] = (a[j] > 0) ? (P[j] -T_r[j])/a[j] : INFINITY;    
    }

    delta = deltat[0];
    mu = 0;
    for (j=0;j<m;j++)
    {
        if (delta > deltat[j])
        {
           delta = deltat[j];
           mu = j;
        }
    }

    t=0;
    for (ti=0;ti<nt;ti++)
    {
        while ((t + delta) < T[ti])
        {
            /*update time*/
            t = t + delta;
            
            /*update state*/
            for (i=0;i<n;i++)
            {
                X[i] += nu[n*mu + i];
            }
            for (j=0;j<m;j++)
            {
                T_r[j] += a[j]*delta;
            }
            P[mu] += durngexps(1);

            /*update propensities*/
            duhzds(m,n,nu_minus,c,X,a);
            for (j=0;j<m;j++)
            {
                deltat[j] = (a[j] > 0 ) ? (P[j] -T_r[j])/a[j] : INFINITY ;    
            }
        
            delta = deltat[0];
            mu = 0;
            for (j=0;j<m;j++)
            {
                if (delta > deltat[j])
                {
                   delta = deltat[j];
                   mu = j;
                }
            }
        }

        /*record our next measurement*/
        for (i=0;i<ndims;i++)
        {
            X_r[i*nt + ti] = X[dims[i]];
        }
    } 

    return 0;
}
