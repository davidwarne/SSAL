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
 * @brief simulate correlated Euler-Maruyama simualtions with nested time
 * time intervals
 * @details Simulates X_f(t) and X_c(t) with t = T[0], ..., T[nt-1] with
 * h_c/h_f = M such that the state space tajectories of X_f and X_c are 
 * strongly correlated.
 *
 * @param m number of parameters
 * @param n dimension of X
 * @param nt number of timesteps
 * @param T array of lenght nt of sampel times
 * @param p vector of model parameters
 * @param X0 the initial condition
 * @param a function pointer for the drift coefficient
 * @param b function pointer for the diffusion coefficient
 * @param ndims number dimension to measure
 * @param dims the dimension indices
 * @param h fine-grain level time-step, the coarse-grain is h*M
 * @param M nesting factor
 * @param X_f_r realisation of X_f(t)
 * @param X_c_r realisation of X_c(t)
 *
 */
int 
dacems(int m,int n, int nt, double * restrict T, double * p, 
    double * restrict X0, 
    void (*a)(double *, unsigned int, double *, unsigned int, double, double *), 
    void (*b)(double*, unsigned int, double *, unsigned int, double, double *), 
    int ndims, int * restrict dims, double h, int M, double * X_f_r, 
    double * X_c_r)
{
    double X_c[n];
    double X_f[n];

    double a_c[n];
    double b_c[n];
    double a_f[n];
    double b_f[n];

    double h_f;
    double h_c;

    double deltaW_c[n];
    double deltaW_f[n];

    double sqrth_f;
    double t_c,t_f;
    

    int i,ti,k;

    h_f = h;
    h_c = h*((double)M);

    /*initial condintions*/
    for (i=0;i<n;i++)
    {
        X_c[i] = X0[i];
    }
    for (i=0;i<n;i++)
    {
        X_f[i] = X0[i];
    }
    for (i=0;i<n;i++)
    {
        deltaW_c[i] = 0;
    }
    for (i=0;i<n;i++)
    {
        deltaW_f[i] = 0;
    }

    t_c = 0;
    t_f = 0;
    k = 0;
    sqrth_f = sqrt(h_f);
    for (ti=0;ti<nt;ti++)
    {
        while ( (t_f + h_f) <= T[ti])
        {
            
            /*compute fine increment*/
            for (i=0;i<n;i++)
            {
                deltaW_f[i] = durngns(0,sqrth_f);
            }
            /*accumulate coarse increment*/
            for (i=0;i<n;i++)
            {
                deltaW_c[i] += deltaW_f[i];
            }
            
            /*compute fine drift and diffusion*/
            (*a)(X_f,n,p,m,t_f,a_f);
            (*b)(X_f,n,p,m,t_f,b_f);

            /* update fine state*/
            for (i=0;i<n;i++)
            {
                X_f[i] += a_f[i]*h_f + b_f[i]*deltaW_f[i];
            }
            t_f += h_f;
            k++;
            
            if (k%M == 0) /*coarse update*/
            {
                /*compute coarse drift and diffusion*/
                (*a)(X_c,n,p,m,t_c,a_c);
                (*b)(X_c,n,p,m,t_c,b_c);

                /* update fine state*/
                for (i=0;i<n;i++)
                {
                    X_c[i] += a_c[i]*h_c + b_c[i]*deltaW_c[i];
                }
                t_c += h_c;
                
                /*reset coarse increment*/
                for (i=0;i<n;i++)
                {
                    deltaW_c[i] = 0;
                }
   
            }
        }

        /*write out data*/
        for (i=0;i<n;i++)
        {
            X_c_r[i*nt + ti] = X_c[dims[i]];
            X_f_r[i*nt + ti] = X_f[dims[i]];
        }
    }
    return 0;    
}
