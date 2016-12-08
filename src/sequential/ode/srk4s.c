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
#include "ODE_sequential.h"

/**
 * @brief 4th order Runge-Kutta method
 * @detail An integration scheme for ordinary differential equations of the form
 *     dY/dt = f(t,Y) where Y \in R^n and f : R^n x [0, \infty) -> R^n
 *
 * @param m the number of parameters,
 * @param n dimension of Y
 * @param nt number of times steps
 * @param T array of nt time points
 * @param p vector of model parameters
 * @param Y0 initial condition
 * @param f vector valued function pointer defining the RHS of the ODE
 * @param ndims number of dimensions to measure
 * @param dims the dimension indices
 * @param h the time step
 * @param Y_r solution stajectory for measured dims (nt*dims)
 */
int srk4s(int m, int n, int nt, float * restrict T, float * restrict p, float * restrict Y0, void (*f)(float *, unsigned int, float *, unsigned int, float, float *), int ndims, int * restrict dims, float h, float * restrict Y_r)
{
    float Y[n]; /*current state*/
    float Y_i[n]; /*intermediate states*/
    float K1[n], K2[n],K3[n], K4[n]; /*RK4 intermediate variables*/
    float t = 0;
    int i, ti;

    for (i=0;i<n;i++)
    {
        Y[i] = Y0[i];
    }

    for (ti=0;ti<nt;ti++)
    {
        for (; t <= (T[ti] -h); t+=h)
        {
            float tj;
            tj = t;
            /*step 1: K1 = f(t,Y)*/
            (*f)(Y,n,p,m,tj,K1);
            /*step 2: K2 = f(t+h/2, Y + K1*h/2)*/
            tj = t + h*0.5; 
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + K1[i]*h*0.5;
            }
            (*f)(Y_i,n,p,m,tj,K2);
            /*step 3: K3 = f(t+h/2, Y + K2*h/2)*/
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + K2[i]*h*0.5;
            }
            (*f)(Y_i,n,p,m,tj,K3);
            /*step 4: K4 = f(t+h, Y + K3*h)*/
            tj = t + h;
            for (i=0;i<n;i++)
            {
                Y_i[i] = Y[i] + K3[i]*h;
            }
            (*f)(Y_i,n,p,m,tj,K4);
            /*complete timestep update*/
            for (i=0;i<n;i++)
            {
                Y[i] += 0.5*h*(K1[i] + K2[i] + K3[i] + K4[i]);
            }
        }
        /*write out timestep*/
        for (i=0;i<ndims;i++)
        {
            Y_r[i*nt+ti] = Y[dims[i]];
        }
    }
    return 0;
}
