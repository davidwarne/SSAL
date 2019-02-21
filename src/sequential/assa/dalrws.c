/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2018  David J. Warne
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
 * @brief Double precision discrete-time, lattice-based random walk with birth
 *        events. 
 * 
 * @param d lattice dimension 
 * @param N array of length d with lattice site counts in each lattice axis. 
 *          size of lattice is \prod_{i=1}^d N[i]
 * @param Nh lattice neighbourhood function (can be used to specify non-square
 *           lattice structures like hex or even irregular lattices).
 * @param Nhmax maximum Neighbourhood size
 * @param nt the number of timesteps to sample at
 * @param T array of length nt of sample times
 * @param p array of model parameters
 * @param X0 initial condition
 * @param g function pointer for crowding function of motility
 * @param f function pointer for crowding function of proliferation
 * @param ndims number of lattice sites to measure
 * @param dims lattice indices to measure
 * @param tau the discrete time-step
 * @param X_r state-space trajectory for measured dims (nt*ndims)
 *
 * @note currently only periodic or reflective boundary conditions can be 
 *        implemented through the Nh function.
 * @note to handle death events or absorbing boundaries etc then the event 
 *       handling must change to allow for population decreases.
 */
int
dalrws(int d, int *N, void (*Nh)(int, int *, int, int *,int *), int Nhmax, int nt, 
       double * restrict T, double * restrict p, double * restrict X0, 
       double (*g)(double, double *), double (*f)(double, double*), 
       int ndims, int * restrict dims, double tau, double * X_r)
{
    int N_tot;
    int X_tot;
    double *X;
    int *Xind;
    double t = 0;
    int *Nhind; 
    
    /* get lattice size */
    N_tot = 1;
    for (int i=0;i<d;i++)
    {
        N_tot *= N[i];
    }

    /*allocate neighbourhood index buffer*/
    Nhind = (int*)malloc(Nhmax*sizeof(int));

    /*allocate buffer memory*/
    X = (double *)malloc(N_tot*sizeof(double));
    Xind = (int *)malloc(N_tot*sizeof(int));

    /*initialise*/
    X_tot = 0;
    for (int i=0;i<N_tot;i++)
    {
        X[i] = X0[i];
        /*build sparse representation*/
        if (X[i] == 1.0)
        {   
            Xind[X_tot++] = i;
        }
    }

    for (int ti=0;ti<nt;ti++)
    {
        for (; t <= (T[ti]-tau);t+=tau)
        {
            int X_tot_prev = X_tot;
            /* perform X_tot trial events*/
            for (int i=0;i<X_tot_prev;i++)
            {
                int Nhsize,j,k;
                double Xnh_avg;
                /*choose a occupied cell at random*/
                j = (int)floor(durngus(0,(double)X_tot));
                /* get the local neighbourhood*/
                (*Nh)(d,N,Xind[j],&Nhsize,Nhind);
                /*choose neighbour random neighbour*/
                k = (int)floor(durngus(0,(double)Nhsize));
                /*attempt event completion if target vacant*/
                if (X[Nhind[k]] == 0.0)
                {
                    double gX,u;
                    /*compute neighbourhood average occupancy*/
                    Xnh_avg = X[Nhind[0]];
                    for (int s=1;s<Nhsize;s++)
                    {
                        Xnh_avg += X[Nhind[s]];
                    }
                    Xnh_avg /= (double)Nhsize;
                    /*compute event completion probabilities*/
                    gX = g(Xnh_avg,p); /*for motility*/
                    u = durngus(0,1);
                    if (u <= gX)
                    {
                        /*move occupancy from j to k*/
                        X[Xind[j]] = 0.0; 
                        X[Nhind[k]] = 1.0;
                        Xind[j] = Nhind[k];
                    }
                }
            }
            /* perform X_tot trial events*/
            for (int i=0;i<X_tot_prev;i++)
            {
                int Nhsize,j;
                double fX,u, Xnh_avg;
                /*choose a occupied cell at random*/
                j = (int)floor(durngus(0,(double)X_tot));
                
                /* get the local neighbourhood*/
                (*Nh)(d,N,Xind[j],&Nhsize,Nhind);
                
                /*compute neighbourhood average occupancy*/
                Xnh_avg = X[Nhind[0]];
                for (int s=1;s<Nhsize;s++)
                {
                    Xnh_avg += X[Nhind[s]];
                }
                Xnh_avg /= (double)Nhsize;
                
                /*compute event completion probabilities*/
                fX = f(Xnh_avg,p); /*for proliferation*/
                u = durngus(0,1);
                if (u <= fX)
                {
                    int k, sum;
                    /*choose unoccupied neighbour random neighbour*/
                    k = (int)floor(durngus(0,(1.0 - Xnh_avg)*((double)Nhsize)));
                    sum = 0;
                    for (int kk=0;kk<Nhsize;kk++)
                    {
                        sum += (X[Nhind[kk]] == 0.0);
                        if (sum == (k + 1))
                        {
                            X[Nhind[kk]] = 1.0;
                            Xind[X_tot++] = Nhind[kk];
                            break;
                        }
                    }
                }
            }
        }

        /* record our next measurement*/
        for (int i=0;i<ndims;i++)
        {
            X_r[i*nt + ti] = X[dims[i]];
        }
    }
    return 0;
}
