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

/**
 * @brief Euler-Maruyama method
 * @detail discretisation scheme for stochastic differential equations 
 * (in the Ito integral form). The method is first order in  weak error
 * and 1/2 order in the strong error. The SDE has the form
 *      dX_t = a(X_t,t)dt + b(X_t,t)dW_t
 * where X_t \in R^n, a() and b() are vector valued functions and W_t is a vector 
 * of independent Wiener processes.
 *
 * The Euler-Maruyama approximation is
 *  X(t + h) = a(X(t),t)h + b(X(t),t)N(0,sqrt(h))
 *
 * @param n dimension of X
 * @param nt numbe of timesteps to sample
 * @param T an array of length nt of sample times
 * @param X0 the initial condition
 * @param a function pointer for the drift coefficient
 * @param b a function pointer for the diffusion coefficient
 * @param ndims numebr dimension to measure
 * @param dims the dimensions indices 
 * @param h the discrete timestep parameter
 * @param X_r state-space trajectory for measured dims (nt*ndims)
 * 
 */
int daems(int n, int nt, double * restrict T, double * restrict X0, void (*a)(double *, int, double,double *), void (*b)(double*,int,double,double*),int ndims, int *restrict dims, double h, double *restrict X_r)
{
    double X[n] /*current state*/
    double ar[n] /* return values of drift function*/
    double br[n] /* return values of the diffusion function*/
    double t=0;
    double deltaW[n]; /*brownian increment*/
    double sqrth;
    
    for (i=0;i<n;i++)
    {
        X[i] = X0[i];
    }
    
    sqrth = sqrt(h);
    for (ti=0;ti<nt;ti++)
    {
        for (;t <= (T[ti] -h ); t+=h)
        {
            /*compute brownian increments*/
            for (i=0;i<n;i++)
            {
                deltaW[i] = durngns(0,sqrth);
            }
            
            /*evaluate drift and diffusion*/
            (*a)(X,n,t,ar);
            (*b)(X,n,t,br);

            /*step in time*/
            for (i=0;i<n;i++)
            {
                X[i] += ar[i]*h + br[i]*deltaW[i];
            }
        }
        /*write out timestep*/
        for (i=0;i<ndims;i++)
        {
            X_r[i*nt+ti] = X[dims[i]];
        }
    }
    return 0;
}
