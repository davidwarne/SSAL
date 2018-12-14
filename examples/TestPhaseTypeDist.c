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
 * @file TestPhaseTypeDist.c
 * @brief Example program for estimating the Phase-type distribution of a 
 * stochastic reaction-diffusion model.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "SSAL.h"

/**
 * @brief Creates 1d instance of reaction-diffusion model in heterogeneous.
 * media.
 * @details Let X be a chemical species and X_i(t) be the copy number of X at 
 * time t in the ith volume of dimension h, i.e., the copy number in 
 * [ih,(i+1)h). Let L = Nh and i = 0,1,...,N-1. 
 *
 *           k1        k1        k2        k2        k2
 *     X_0 <---> X_1 <---> X_2 <---> ... <---> X_N-1 --> A
 * |<---h--->|<---h--->|<---h--->|   ...   |<---h--->)
 * |<--------------------->|<----------------------->|
 * 0                       S                         L
 * where k1 = D1/h^2 where D1 is the diffusivity in the continuum limit in [0,S)
   k2 = D2/h^2 where D2 is the diffusivity in the coninuum limit in [S,L) and A is 
 * an absorptive boundary at L. X_0 implicitly includes a reflective boundary.
 *
 *
 * @todo would like to add a Reaction Diffusion structure to simplify this 
 * long term a bit. could have: 1) a CRN; 2) diffusivity function;
 * and 3) boundaries/interfaces. The some functions to convert to compartment-based
 * RDME or SDE/indiviual based model.
 *
 * @todo extend to two regions with interface cell location
 */
SSAL_CRN 
CreateSpatialCRN(SSAL_real_t L1, SSAL_real_t L2, SSAL_real_t h, SSAL_real_t D1, SSAL_real_t D2)
{
    int N,N1; /*number of cells*/
    SSAL_CRN newCRN;
    char ** names;          /*strings of printable variable names*/
    SSAL_real_t * nu_minus; /* nu^- reactant stoichiometric coefficients*/
    SSAL_real_t * nu_plus;  /* nu^+ product stoichiomeric coefficients*/
    SSAL_real_t * c;        /* rate parameters c[i] = k for all i here*/

    N = (int)ceil(L2/h); /* adjust h if needed to ensure L = Nh and N is int*/
    if (L2 < ((SSAL_real_t)N)*h)
    {
        h = L2/((SSAL_real_t)N);
    }

    /* build arrays required for the CRN*/
    names = (char**)malloc((N+1)*sizeof(char*));
    names[0] = (char*)malloc((N+1)*SSAL_MAX_NAME_SIZE*sizeof(char));
    for (int i=1;i<N+1;i++)
    {
        names[i] = names[i-1] + SSAL_MAX_NAME_SIZE;
    }
    nu_minus = (SSAL_real_t *) malloc((N+1)*(2*N -1)*sizeof(SSAL_real_t));
    nu_plus = (SSAL_real_t *) malloc((N+1)*(2*N -1)*sizeof(SSAL_real_t));
    c = (SSAL_real_t *) malloc((2*N -1)*sizeof(SSAL_real_t));
    
    /* set names */
    for (int i=0;i<N;i++)
    {
        snprintf(names[i],SSAL_MAX_NAME_SIZE,"X_%d",i+1); 
    }
    snprintf(names[N],SSAL_MAX_NAME_SIZE,"A");
    
    /* set stoichiometries*/
    memset(nu_minus,0,(N+1)*(2*N+1)*sizeof(SSAL_real_t));
    memset(nu_plus,0,(N+1)*(2*N+1)*sizeof(SSAL_real_t));
    nu_minus[0] = 1.0;
    nu_plus[1] = 1.0;
    for (int i=1;i<N;i++)
    {
        /* set reactant stoichiometries*/
        nu_minus[(2*i-1)*(N+1) + i] = 1.0; 
        nu_minus[(2*i)*(N+1) + i] = 1.0;
        /* set product stoichiometries*/
        nu_plus[(2*i-1)*(N+1) + i -1] = 1.0; 
        nu_plus[(2*i)*(N+1) + i + 1] = 1.0;
    }
    /*rate parameters*/
    /*find the cell that contains the interface S*/
    N1 = (int)floor(L1/h);
    for (int i=0;i<2*N1-1;i++)
    {
        /* compute rate constant to ensure diffusivity D1 is recovered for h -> 0*/
        c[i] = D1/(h*h);
    }
    for (int i=2*N1;i<2*N-1;i++)
    {
        /* compute rate constant to ensure diffusivity D1 is recovered for h -> 0*/
        c[i] = D2/(h*h);
    }
    /* construct the CRN*/
    newCRN = SSAL_CreateChemicalReactionNetwork(names,2*N-1,N+1,nu_minus,nu_plus,c);
    free(names[0]);
    free(names);
    free(nu_minus);
    free(nu_plus);
    free(c);
    return newCRN;
}


int
main(int argc, char ** argv)
{
    /* spatial model paramaters*/
    double D[2],L[2],h,dt;
    int NT, NR;
    SSAL_real_t *T;
    SSAL_real_t *A; /*only record the absorptive boundary*/
    SSAL_CRN CRN;
    /* get input args */
    for (int i=0;i<argc;i++)
    {
        if (!strcmp("-n", argv[i]))
        {
            NR = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-p",argv[i]))
        {
            D[0] = (double)atof(argv[++i]);
            D[1] = (double)atof(argv[++i]);
            L[0] = (double)atof(argv[++i]);
            L[1] = (double)atof(argv[++i]);
            h = (double)atof(argv[++i]);
        }
        else if (!strcmp("-nt",argv[i]))
        {
            NT = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-t",argv[i]))
        {
            if (NT == 0)
            {
                fprintf(stderr,"Must use the -nt option to choose number of time points\n");
                exit(1);
            }
            T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
            T[NT-1] = (SSAL_real_t)atof(argv[++i]);
            /*build sample times*/
            dt = T[NT-1]/((SSAL_real_t)(NT));
            for (i=0; i < (NT-1); i++){
                T[i] = dt*((SSAL_real_t)(i+1));
            }
        }
        else if (!strcmp("-h",argv[i]))
        {
            fprintf(stderr,"Usage: %s -n numRealisations -p D1 D2 L1 L2 h -nt numTimePoints -t endTime",argv[0]);
            exit(0);
        }
    }
    
    /* init SSAL*/
    SSAL_Initialise(argc,argv);

    /* create the spatial CRN*/
    CRN = CreateSpatialCRN(L[0],L[1], h, D[0],D[1]);
    /* since NR single particle realisation is equivalent to 1 realisation with X_0(0) = NR*/
    CRN.X0 = (SSAL_real_t *)malloc(CRN.N*sizeof(SSAL_real_t));
    memset(CRN.X0,0,CRN.N*sizeof(SSAL_real_t));
    CRN.X0[0] = (SSAL_real_t)NR;
    /*only recording  the absorptive boundary*/
    CRN.nvar = 1;
    CRN.vars = (int *)malloc(CRN.nvar*sizeof(int));
    CRN.vars[0] = CRN.N-1;

    A = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
    /*run simulation*/
    demnrms(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,CRN.nu,CRN.c,CRN.nvar,CRN.vars,A);

    fprintf(stdout,"\"t\",\"A(t)\"\n");
    for (int i=0;i<NT;i++)
    {
        fprintf(stdout,"%g,%g\n",T[i],A[i]);
    }


    exit(0);
}
