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
 * @brief simulate correlated tau-leaping simulations with nested times intervals
 * @details Simulates Z_l(t) and Z_{l-1}(t) with t = T[0], ... , T[nt-1] with 
 * tau_{l-1}/tau_l = M such that the state space trajectories of Z_l(t) 
 * and Z_{l-1}(t) are strongly correlated.
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
 * @param tau fine-grain level time-step coarse-grain is tau*M;
 * @param M nesting factor
 * @param Z_l_r realisation of Z_l(t)
 * @param Z_lm1_r realisation of Z_{l-1}(t) with strong correlation to Z_l(t)
 */
int 
dactauls(int m, int n, int nt, double * restrict T, double * restrict X0, 
         double * restrict nu_minus, double * restrict nu, double * restrict c, 
         int ndims , int * restrict dims,double tau, int M, 
         double * restrict Z_f_r, double * restrict Z_c_r)
{
    double Z_c[n]; /*coarse and fine grain state vectors*/ 
    double Z_f[n];
    double a_c[m]; /*coarse and fine grain propensity values*/
    double a_f[m]; 

    /*virtual propensity channels*/
    double b[3][m];
    /*poisson variates*/
    double Y[3][m];
    double t = 0;
    int ti,i,j,k,r;
   
    double tau_f;
    double tau_c;
    double t_c;
    double t_f;

    tau_f = tau;
    tau_c = tau*((double)M);
    
    /*initial conditions*/
    for (i=0;i<n;i++)
    {
        Z_c[i] = X0[i];
    }

    for (i=0;i<n;i++)
    {
        Z_f[i] = X0[i];
    }

    
    /*compute coarse propensities*/
    duhzds(m,n,nu_minus,c,Z_c,a_c);
    /*compute fine propensities*/
    duhzds(m,n,nu_minus,c,Z_f,a_f);

    t_c = 0;
    t_f = 0;

    for (ti=0;ti<nt;ti++)
    {
        while((t_f + tau_f) <= T[ti])
        {
            /*compute fine propensities*/
            duhzds(m,n,nu_minus,c,Z_f,a_f);
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
                    /*we only generate a Poisson RV for b > 0*/
                    Y[r][j] = (b[r][j] <= 0) 
                              ? 0 : (double)durngpois(b[r][j]*tau_f);
                }
            }

            /*update coarse and fine state vectors*/
            /*note po(a_c*tau_l) = po(b_1*tau_l) + po(b_2*tau_l)*/
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
            t_f += tau_f;
            k++;

            if (k%M == 0) /*coarse update*/
            {
                /*compute coarse propensities*/
                duhzds(m,n,nu_minus,c,Z_c,a_c);
                t_c += tau_c;
            }

        }
        /*write out timesteps*/
        for (i=0;i<ndims;i++)
        {
            Z_f_r[i*nt+ti] = Z_f[dims[i]];
        }
        for (i=0;i<ndims;i++)
        {
            Z_c_r[i*nt+ti] = Z_c[dims[i]];
        }
    }

    return 0;
}
