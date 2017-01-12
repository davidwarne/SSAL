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
 * @file TestSSAL.c
 * @brief An example program using the SSAL API to sample a discrete
 * state continuous time markov process
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include<stdlib.h>
#include<stdio.h>
#include "SSAL.h"

#define EXACT_COUPLING 0
#define COUPLED_TAU_LEAP 1

int 
main(int argc , char ** argv)
{
    int i,j,k;
    int NT;
    int NR; /*number of realisations*/
    SSAL_real_t dt;
    SSAL_real_t *T;
    SSAL_real_t *X0;
    SSAL_real_t * X_f_r;
    SSAL_real_t * X_c_r;
    SSAL_real_t * nu;
    char **names;
    int type;
    char opts[256];
    int m,n;
    int M;
    /*tau-leap args*/
    SSAL_real_t tau;
    
    char *filename;
    int algType;
    SSAL_CRN CRN;
    int *d; 
    NT = 0;    
    NR = 1;
    algType = EXACT_COUPLING;
    opts[0] = '\0';

    /*get input args with getops*/
    for (i=0;i<argc;i++)
    {
        if (!strcmp("-n",argv[i]))
        {
            NR = (int)atoi(argv[++i]);
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
        }
        else if (!strcmp("-f",argv[i]))
        {
            filename = argv[++i];
        }
        else if (!strcmp("-st",argv[i]))
        {
            type = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-T",argv[i]))
        {
            tau = (SSAL_real_t)atof(argv[++i]);
        }
        else if (!strcmp("-M",argv[i]))
        {
            algType = COUPLED_TAU_LEAP;
            M = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-h",argv[i]))
        {
            fprintf(stderr,"Usage: %s [-n numRealisations] -nt numTimePoints -t endTime -f filename [-st simType] [-T tau] [-M M]\n",argv[0]);
            exit(0);
        }
    }
    
    /*build sample times*/
    dt = T[NT-1]/((SSAL_real_t)(NT));
    for (i=0; i < (NT-1); i++){
        T[i] = dt*((SSAL_real_t)(i+1));
    }

    

    /*init SSAL */
    SSAL_Initialise(argc,argv);
    /*import model*/
    CRN = SSAL_ImportLSBML(filename);
    /*Sanity check*/ 
    SSAL_WriteChemicalReactionNetwork(stderr,CRN);
    
    /*we will just use the default initial conditions*/
    X0 = CRN.X0;
    names = CRN.names;
    n = CRN.N;
    for(i=0;i<CRN.N;i++)
        fprintf(stderr,"%s\n",names[i]);

    d = (int *)malloc(CRN.N*sizeof(int));
    nu = (SSAL_real_t *)malloc(CRN.N*CRN.M*sizeof(SSAL_real_t));
    X_f_r = (SSAL_real_t *)malloc(CRN.N*NT*NR*sizeof(SSAL_real_t));
    X_c_r = (SSAL_real_t *)malloc(CRN.N*NT*NR*sizeof(SSAL_real_t));

    for (i=0;i<CRN.N*CRN.M;i++)
    {
        nu[i] = CRN.nu_plus[i] - CRN.nu_minus[i];
    }
    for (i=0;i<CRN.N;i++)
    {
        d[i] = i; 
    }

    switch (algType)
    {
        default:
        case EXACT_COUPLING:
            for(i=0;i<NR;i++)
            {
                daectauls(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,nu,CRN.c,CRN.N,d,
                       tau,X_f_r+i*CRN.N*NT,X_c_r + i*CRN.N*NT);
            }
            break;
        case COUPLED_TAU_LEAP:
            for(i=0;i<NR;i++)
            {
                dactauls(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,nu,CRN.c,CRN.N,d,
                       tau,M,X_f_r+i*CRN.N*NT,X_c_r+i*CRN.N*NT);
            }

            break;
    }

    /*write the data*/
    fprintf(stdout,"0");
    for (j=0;j<NT;j++)
    {
        fprintf(stdout,",%g",T[j]);
    }
    fprintf(stdout,"\n");
    for (i=0;i<NR;i++)
    {
        for (k=0;k<CRN.N;k++)
        {
            fprintf(stdout,"%f",CRN.X0[k]);
            for (j=0;j<NT;j++)
            {
                fprintf(stdout,",%g",X_f_r[(i*CRN.N + k)*NT +j]);
            }
            fprintf(stdout,"\n");
        }
        for (k=0;k<CRN.N;k++)
        {
            fprintf(stdout,"%f",CRN.X0[k]);
            for (j=0;j<NT;j++)
            {
                fprintf(stdout,",%g",X_c_r[(i*CRN.N + k)*NT +j]);
            }
            fprintf(stdout,"\n");
        }

    }
    return 0;
}
