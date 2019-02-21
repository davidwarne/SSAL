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
/**
 * @file TestDRW.c
 * @brief example of a discrete-time lattice based random walk model.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 */
#include <stdlib.h>
#include <stdio.h>
#include "SSAL.h"

/* the hexagonal lattice model of Jin et al. 2016*/

/* neighbourhood function for the lattice*/
void
hexNh(int d, int *N, int ind, int *Nhsize, int *hex)
{
    /* we assume d = 2, this paramter is intendeded for future generalisations*/
    int Nx = N[0];
    int Ny = N[1];
    int offx[6] = {0,-1,0,1,1,1};
    int offy[6] = {-1,0,1,1,0,-1};
    /* get 2d index*/
    int j = ind / Nx;
    int i = ind % Nx;

    *Nhsize = 6;
    /*implementing periodic boundary*/
    /* note, for periodic boundary to be valid Ny must be even*/
    for (int k=0;k<*Nhsize;k++)
    {   
        hex[k] = ((j+offy[k]+Ny)%Ny)*Nx + (i+offx[k]+Nx)%Nx;
    }
}

/* motility function*/
double g(double Cavg, double *params)
{
    double Pm = params[0];
    return Pm; /*leads to linear diffusion in continuum limit*/
}

/* proliferation function*/
double f(double Cavg, double *params)
{
    double Pp = params[1];
    double alpha = params[2];
    double beta = params[3];
    double gamma = params[4];
    /* leads to generalised logistic growth in the continuum limit 
     * dC/dt = \lambda C^\alpha (1-C^\beta)^\gamma
     */
    //fprintf(stderr,"f(c) =  %g\n",Pp*pow(Cavg,alpha-1.0)*pow(1.0 - pow(Cavg, beta),gamma)/(1.0-Cavg));
    return Pp*pow(Cavg,alpha-1.0)*pow(1.0 - pow(Cavg, beta),gamma); /*leads to logistic growth in continuum limit*/
    //return Pp; /*leads to logistic growth in continuum limit*/
    
}

int main(int argc, char ** argv)
{
    SSAL_real_t tau;
    SSAL_real_t Pm;
    SSAL_real_t Pp;
    SSAL_real_t alpha;
    SSAL_real_t beta;
    SSAL_real_t gamma;
    SSAL_real_t *C0;
    SSAL_real_t *T;
    SSAL_real_t params[5];
    SSAL_real_t* C_r;
    SSAL_real_t* C;
    int N[2];
    int nt;
    int I,J;
    int *dims;
    int NR;

    /*fixed parameters*/
    Pm = 1; tau = 1.0;
    if (argc < 9)
    {
        fprintf(stderr,"Usage %s [I] [J] [timesteps] [NR] [Pp] [alpha] [beta] [gamma]\n",argv[0]);
        exit(1);
    }
    /*get time steps and Proliferation parameters*/
    I = (int)atoi(argv[1]);
    J = (int)atoi(argv[2]);
    nt = (int)atoi(argv[3]);
    NR = (int)atoi(argv[4]);
    Pp = (SSAL_real_t)atof(argv[5]);
    alpha = (SSAL_real_t)atof(argv[6]);
    beta = (SSAL_real_t)atof(argv[7]);
    gamma = (SSAL_real_t)atof(argv[8]);


    N[0] = I; N[1] = J;
    params[0] = Pm; params[1] = Pp;
    params[2] = alpha; params[3] = beta; params[4] = gamma;

    T = (SSAL_real_t*)malloc(nt*sizeof(SSAL_real_t));
    for (int j=0;j<nt;j++)
    {
        T[j] = (SSAL_real_t)(j+1);
    }

    /* init the RNG */
    SSAL_Initialise(argc,argv);
    
    /*initialise lattice and initial occupacies*/
    C0 = (SSAL_real_t*) malloc(I*J*sizeof(SSAL_real_t));
    dims = (int*)malloc(I*J*sizeof(int));
    for (int i =0;i<I*J;i++)
    {
        SSAL_real_t u;
        u = durngus(0,1);
        C0[i] =  (SSAL_real_t)(u <= 0.05);
        dims[i] = i; /*observe all sites*/
    }

    C_r = (SSAL_real_t*)malloc(I*J*nt*sizeof(SSAL_real_t));
    memset((void*)C_r,0,I*J*nt*sizeof(SSAL_real_t));
    C = (SSAL_real_t*)malloc(nt*sizeof(SSAL_real_t));
    memset((void*)C,0,nt*sizeof(SSAL_real_t));

    /*simulate realisations*/
    for (int i=0;i<NR;i++)
    {
        dalrws(2,N,&hexNh,6,nt,T,params,C0,&g,&f,I*J,dims,tau,C_r);
        for (int i=0;i<I*J;i++)
        {
            for (int ti=0;ti<nt;ti++)
            {
                C[ti] += C_r[i*nt + ti];        
            }
        }
    }

    /*take spatial average*/
    for (int ti=0;ti<nt;ti++)
    {
        C[ti] /= (SSAL_real_t)(I*J*NR);        
    }

    fprintf(stdout,"\"t\",\"C(t)\"\n");
    for (int ti=0;ti<nt;ti++)
    {
        fprintf(stdout,"%d,%f\n",ti,C[ti]);
    }
    free(C0);
    free(C_r);
    free(C);
    free(dims);
    exit(0);

}
