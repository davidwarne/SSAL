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

#include "ASSA_sequential.h"
#include "util_sequential.h"

/**
 * @brief Simulated correlated tau-leap and exact samples paths.
 * @details Generates Z_tau(t) and X(t) with t = T[0],...,T[nt-1] 
 * with Z_tau being an Euler approximation of X(t) such that the 
 * state space trajectories Z_tau(t) and X(t) are strongly correlated.
 *
 * @param m the number of reactions
 * @param n the dimension of the state space
 * @param nt number of timesteps to sample at
 * @param T the sample times
 * @param X0 initial conditions
 * @param nu_minus reactant stoichiometric coefficients m x n (row-major format)
 * @param nu the stoichiometric matrix m x n (row-major format)
 * @param c kinetic reaction rate parameters
 * @param ndims number of dimensions to measure
 * @param dims vector of dimension indices indicating those to be measured
 * @param tau discretisation step for tau-leaping meathod
 * @param Z_r realisation of Z(t)
 * @param X_r realisation of X(t)
 *
 * @note the coupling is based on Anderson and Higham's method which couples a
 * tau-leaping simulation with a modified next-reaction simulation.
 */
int daectauls(int m, int n, int nt, double *restrict T, double * restrict X0,
    double * restrict nu_minus, double *restrict nu, double * restrict c, int ndims,
    int * dims, double tau, int M, double * restrict Z_r, double * restrict X_r)
{
    double Z[n]; /*approximate and exact states*/
    double X[n];
    double Ztilde[n]; /*stored approximate state of last step*/
    double a_Z[n]; /*exact and approximate propensities*/
    double a_X[n];
    /*virual propensity channels*/
    double b[3][m];
    /*P ~ Exp(1) */
    double P[3][m];
    double T_r[3][m];
    double deltat[3][m];
    int mu[2];

    int i,j,k,r,ti;
    double t, delta, T_tau;

    /*initialise*/
    for (i=0;i<n;i++)
    {
        X[i] = X0[i];
    }
    for (i=0;i<n;i++)
    {
        Z[i] = X0[i];
    }
    
    t = 0;
    
    for (i=0;i<n;i++)
    {
        Ztilde[i] = Z[i];
    }

    T_tau = tau;
    for (r=0;r<3;r++)
    {
        for (j=0;j<m;j++)
        {
            P[r][j] = durngexps(1);
        }
    }
    for (r=0;r<3;r++)
    {
        for (j=0;j<m;j++)
        {
            T_r[r][j] = 0;
        }
    }

    /*compute propensities*/
    duhzds(m,n,nu_minus,c,Ztilde,a_Z);
    duhzds(m,n,nu_minus,c,X,a_X);
   
    /*both propensities will be the same at t = 0*/
    for (j=0;j<m;j++)
    {
        b[0][j] = a_X[j];
        b[1][j] = 0;
        b[2][j] = 0;
    }
   
    for (r=0;r<3;r++)
    {
        for (j=0;j<m;j++)
        {
            deltat[r][j] = (b[r][j] > 0) ? (P[r][j] - T_r[r][j])/b[r][j] : INFINITY;
        }
    }
    /*only need to scan r = 0*/
    delta = delta[0][0];
    for (j=0;j<m;j++)
    {
        if (delta > delta[0][j])
        {
            delta = delta[0][j];
            mu[0] = 0;
            mu[1] = j;
        }
    }

    for (ti=0;ti<nt;ti++)
    {
        while((t+delta) <= T[ti])
        {
            if (t + delta >= T_tau)
            {
                /*store current state for propensity update*/
                for (i=0;i<n,i++)
                {
                    Ztilde[i] = Z[i];
                }
                for (r=0;r<3;r++)
                {
                    for (j=0;j<m;j++)
                    {
                        T_r[r][j] += b[r][j]*(T_tau - t);
                    }
                }
                t = T_tau;
                T_tau += tau;
            }
            else
            {
                /*update states*/
                for(i=0;i<n;i++)
                {
                    X[i] += (mu[0] != 2) ? nu[n*mu[1] +i] : 0;
                }
                for (i=0;i<n;i++)
                {
                    Z[i] += (mu[0] != 1) ? nu[n*mu[1] +i] : 0;
                }
                
                for (r=0;r<3;r++)
                {
                    for (j=0;j<m;j++)
                    {
                        T_r[r][j] += b[r][j]*delta;
                    }
                }

                P[mu[0]][mu[1]] += durngexps(1);
                t = t + delta;
            }

            /* update propensities*/
            duhzds(m,n,nu_minus,c,Ztilde,a_Z);
            duhzds(m,n,nu_minus,c,X,a_X);
            
            /*update virtual propensities*/
            for (j=0;j<m;j++)
            {
                b[0][j] = (a_Z[j] < a_X[j]) ? a_Z[j] : a_X[j];
            }
            for (j=0;j<m;j++)
            {
                b[1][j] = a_X[j] - b[0][j];
            }
            for (j=0;j<m;j++)
            {
                b[2][j] = a_Z[j] - b[0][j];
            }
            /*update deltas*/
            for (r=0;r<3;r++)
            {
                for (j=0;j<m;j++)
                {
                    deltat[r][j] = (b[r][j] > 0) ? (P[r][j] - T_r[r][j])/b[r][j] : INFINITY;
                }
            }
            /*find next reaction*/
            delta = deltat[0][0];
            for (r=0;r<3;r++)
            {
                for (j=0;j<m;j++)
                {
                    if (delta > deltat[r][j])
                    {
                        delta = delta[r][j];
                        mu[0] = r;
                        mu[1] = j;
                    }
                }
            }
        }
        /*write out timesteps*/
        for (i=0;i<ndims;i++)
        {
            X_r[i*nt + ti] = X[dims[i]];
        }
        for (i=0;i<ndims;i++)
        {
            Z_r[i*nt + ti] = Z[dims[i]];
        }
    }
}
