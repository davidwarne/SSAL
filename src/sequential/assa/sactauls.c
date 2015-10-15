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
 * @brief simulate correlated tau-leaping simulations with nested times intervals
 * @details Simulates Z_l(t) and Z_{l-1}(t) with t = T[0], ... , T[nt-1] with 
 * tau_{l-1}/tau_l = M such that the state space trajectories of Z_l(t) and Z_{l-1}(t)
 * are strongly correlated.
 *
 * @param m the number of reactions
 * @param n the dimension of state space
 * @param nt number of timesteps to sample at
 * @param T the sample times
 * @param X0 initial conditions
 * @param nu_minus reactant stoichiometric coefficients m x n (row-major format)
 * @param nu stoichiometric matrix m x n (row-major format)
 * @param c kinetic reaction rates
 * @param ndims number of dimensions to measure 
 * @param dims vector of dimesnion indices indicating those to be mesured
 * @param tau coarse-grain level time-step fine-grain is tau/M;
 * @param M nesting factor
 * @param Z_l_r realisation of Z_l(t)
 * @param Z_lm1_r realisation of Z_{l-1}(t) with strong correlation to Z_l(t)
 */
int sactauls(int m, int n, int nt, float * restrict T, float * restrict X0, 
   float * restrict nu_minus, float * restrict nu, float * restrict c, int ndims , 
   int * restrict dims,float tau, int M, float * restrict Z_l_r, float * restrict Z_lm1_r)
{
    float Z_c[n]; /*coarse and fine grain state vectors*/ 
    float Z_f[n];
    float a_c[m]; /*coarse and fine grain propensity values*/
    float a_f[m]; 

    /*virtual propensity channels*/
    float b[3][m];
    /*poisson variates*/
    float Y[3][m];
    float t = 0;
    int ti,i,j,k,r;
   
    float tau_l;
    float tau_lm1;

    tau_l = tau;
    tau_lm1 = tau*M;
    /*initial conditions*/
    for (i=0;i<n;i++)
    {
        Z_c[i] = X0[i];
    }

    for (i=0;i<n;i++)
    {
        Z_f[i] = X0[i];
    }


    for (ti=0;ti<nt;ti++)
    {
        for (;t<=T[ti];t+=tau_lm1)
        {
            /*compute coarse propensities*/
            suhzds(m,n,nu_minus,c,Z_c,a_c);
            
            for (k=0;k<M;k++)
            {
                /*compute fine propensities*/
                suhzds(m,n,nu_minus,c,Z_f,a_f);
                /*update virtual propensities*/
                for (j=0;j<m;j++)
                {
                    b[0][j] = (a_c[j] < a_f[j]) ? a_c[j] : a_f[j];
                }
                for (j=0;j<m;j++)
                {
                    b[1][j] = a_c[j] - b[0][j];
                }
                for (j=0;j<m;j++)
                {
                    b[2][j] = a_f[j] - b[0][j];
                }
                /*generate poisson variates for each virtual reaction channel*/
                for (r=0;r<3;r++)
                {
                    for (j=0;j<m;j++)
                    {
                        Y[r][j] = (b[r][j] <= 0) ? 0 : (float)surngpois(b[r][j]*tau_l);
                    }
                }


                /*update coarse and fine state vectors*/
                /*note Po(a_c*tau_l) = Po(b_1*tau_l) + Po(b_2*tau_l)*/
                for (j=0;j<m;j++)
                {
                    for (i=0;i<n;i++)
                    {
                        Z_c[i] += (Y[0][j] + Y[1][j])*nu[j*n+i];
                    }
                }
                
                /*note Po(a_f*tau_l) = Po(b_1*tau_l) + Po(b_3*tau_l)*/
                for (j=0;j<m;j++)
                {
                    for (i=0;i<n;i++)
                    {
                        Z_f[i] += (Y[0][j] + Y[2][j])*nu[j*n+i];
                    }
                }

            }
        }
        /*write out timesteps*/
        for (i=0;i<ndims;i++)
        {
            Z_l_r[i*nt+ti] = Z_f[dims[i]];
        }
        for (i=0;i<ndims;i++)
        {
            Z_lm1_r[i*nt+ti] = Z_c[dims[i]];
        }
    }

    return 0;
}
