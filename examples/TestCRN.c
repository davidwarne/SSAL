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

#define GILLESPIE 0
#define MNRM 2
#define TAU_LEAP 1

int 
main(int argc , char ** argv)
{
    int i,j,k;
    int NT;
    int NR; /*number of realisations*/
    SSAL_real_t dt;
    SSAL_real_t *T;
    SSAL_real_t *X0;
    SSAL_real_t * X_r;
    char **names;
    int type;
    char opts[256];
    int m,n;

    /*tau-leap args*/
    SSAL_real_t tau;
    
    char *filename;
    int algType;
    SSAL_CRN CRN;
    NT = 0;    
    NR = 1;
    algType = MNRM;
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
            algType = TAU_LEAP;
            tau = (SSAL_real_t)atof(argv[++i]);
        }
        else if (!strcmp("-h",argv[i]))
        {
            fprintf(stderr,"Usage: %s [-n numRealisations] -nt numTimePoints -t endTime -f filename [-st simType] [-T tau] [-T0 tau0 -L L -M M -eps eps]\n",argv[0]);
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

    X_r = (SSAL_real_t *)malloc(CRN.N*NT*NR*sizeof(SSAL_real_t));

    switch (algType)
    {
        default:
        case MNRM:
            for(i=0;i<NR;i++)
            {
            //    degils(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,CRN.nu,CRN.c,
            //           CRN.nvar,CRN.vars,X_r+i*CRN.N*NT);
                demnrms(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,CRN.nu,CRN.c,
                       CRN.nvar,CRN.vars,X_r+i*CRN.N*NT);
            }
            break;
        case TAU_LEAP:
            for(i=0;i<NR;i++)
            {
                datauls(CRN.M,CRN.N,NT,T,CRN.X0,CRN.nu_minus,CRN.nu,CRN.c,
                        CRN.nvar,CRN.vars,tau,X_r+i*CRN.N*NT);
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
                fprintf(stdout,",%g",X_r[(i*CRN.N + k)*NT +j]);
            }
            fprintf(stdout,"\n");
        }
    }
    return 0;
}
