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
/**
 * @file TestPDE.c
 * @brief An example program using the PDE solver
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "SSAL.h"

#define K_CAR_CAP 1.7e-3

/*example of PDEs using Non-linear BTCS*/

/*Fisher-KPP equation*/

/* dc/dt = Dd^2c/dc^2 + \lambda c(1-c)*/
void
D_FKPP(double *c, int n, double *params, double *D)
{
    /*linear diffusion, so just a constant*/
    double D_0;
    D_0 = params[0];
    for (int i=0;i<n;i++)
    {
        D[i] = D_0;
    }
}

/*Porous-Fisher equation*/
/* dc/dt = D_0 dJ/dx + \lambda c(1-c) with J = c dc/xu*/
void
D_PF(double *c, int n, double *params, double *D)
{
    double D_0;
    D_0 = params[0];
    for (int i=0;i<n;i++)
    {
        D[i] = D_0*c[i]/K_CAR_CAP;
    }
}

/*logistic growth*/
void
lgf(double *c, int n,double *params, double *S)
{
    double lambda;
    lambda = params[1];
    for (int i=0;i<n;i++)
    {
        S[i] = lambda*c[i]*(1.0 - c[i]/K_CAR_CAP);
    }
}

int 
main(int argc , char ** argv)
{
    double T;
    double L;
    double *Ti;
    int NX;
    int NT;
    double tol;
    double D_0;
    double lambda;
    double *C0;
    double *C_r;
    double *params;
    clock_t start_t,end_t;

    // crazy large parameters to test functionality
    D_0 = 5000.0;
    lambda = 1;
    NT = 16000;
    tol = 1e-6;
    
    if (argc < 5)
    {
        fprintf(stderr,"Usage: %s D lambda NT tol\n",argv[0]);
        exit(1);
    }
    else
    {
        D_0 = (double)atof(argv[1]);
        lambda = (double)atof(argv[2]);
        NT = (int)atoi(argv[3]);
        tol = (double)atof(argv[4]);
    }

    SSAL_Initialise(argc,argv);
    
    //D_0 = 1e-3;
    //lambda = 4e-3;
    L = 1875.0;
    T = 48.0;
    NX = 186;

    C0 = (double*)malloc(NX*sizeof(double));
    memset(C0,0,NX*sizeof(double));
    for (int i=0;i<20;i++)
    {
        C0[i] = 0.7*K_CAR_CAP;
    }
   
    for (int i=156;i<186;i++)
    {
        C0[i] = 0.7*K_CAR_CAP;
    }

    Ti = (double*)malloc(NT*sizeof(double));
    for(int i=0;i<NT;i++)
    {
        Ti[i] = ((double)i)*T/((double)(NT-1));
    }

    C_r = (double *)malloc(NT*NX*sizeof(double));
    memset(C_r,0,NX*NT*sizeof(double));

    params= (double*)malloc(2*sizeof(double));
    params[0] = D_0;
    params[1] = lambda;
    
    start_t = clock();
    if(dbtcsfps(L,T,NT,Ti,params,C0,&D_PF,&lgf,NX,NT,tol,1000,C_r) < 0)
    {
        free(C0);
        free(C_r);
        exit(1);
    }
    end_t = clock();

    fprintf(stderr,"time = %f seconds\n", ((double)(end_t-start_t))/((double)CLOCKS_PER_SEC));

    for (int j=0;j<NT;j++)
    {
        fprintf(stdout,"%f",Ti[j]);
        for (int i=0;i<NX;i++)
        {
            fprintf(stdout,",%f",C_r[j*NX + i]);
        }
        fprintf(stdout,"\n");
    }
    free(C0);
    free(C_r);
    exit(0); 
}
