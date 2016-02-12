/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2015  David J. Warne
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
 * @date 20 Sep 2015
 */
#include<stdlib.h>
#include<stdio.h>
#include "SSAL.h"


int main(int argc , char ** argv)
{
    int i;
    int NT;
    int NR; /*number of realisations*/
    float dt;
    float *T;
    float *X0;
    char **names;
    //int model;
    int type;
    char opts[256];
    int m,n;

    /*tau-leap args*/
    float tau;
    /*Multi-level args*/
    int M,L;
    float tau0,eps;

    char *filename;
    SSAL_AlgorithmType algType;
    SSAL_Model CRN;
    SSAL_ChemicalReactionNetwork *CRN_ptr;
    SSAL_Simulation sim;
    type = 0;
    NT = 0;    
    NR = 1;
    algType = SSAL_ESSA_GILLESPIE_SEQUENTIAL;
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
            T = (float *)malloc(NT*sizeof(float));
            T[NT-1] = (float)atof(argv[++i]);
        }
        else if (!strcmp("-f",argv[i]))
        {
            filename = argv[++i];
        }
//        else if (!strcmp("-m",argv[i]))
//        {
//            model = (int)atoi(argv[++i]);
//        }
        else if (!strcmp("-st",argv[i]))
        {
            type = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-T",argv[i]))
        {
            algType = SSAL_ASSA_TAU_LEAP_SEQUENTIAL;
            tau = (float)atof(argv[++i]);
        }
        else if(!strcmp("-L",argv[i]))
        {
            algType = SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL;
            L = (int)atoi(argv[++i]);
        }
        else if(!strcmp("-M",argv[i]))
        {
            algType = SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL;
            M = (int)atoi(argv[++i]);
        }
        else if(!strcmp("-T0",argv[i]))
        {
            algType = SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL;
            tau0 = (float)atof(argv[++i]);
        }
        else if(!strcmp("-eps",argv[i]))
        {
            algType = SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL;
            eps = (float)atof(argv[++i]);
        }
        else if (!strcmp("-h",argv[i]))
        {
            //fprintf(stderr,"Usage: %s [-n numRealisations] -nt numTimePoints -t endTime [-m model] [-st simType] [-T tau] [-T0 tau0 -L L -M M -eps eps]\n");
            fprintf(stderr,"Usage: %s [-n numRealisations] -nt numTimePoints -t endTime -f filename [-st simType] [-T tau] [-T0 tau0 -L L -M M -eps eps]\n",argv[0]);
            exit(0);
        }
    }
    
    /*build sample times*/
    dt = T[NT-1]/((float)(NT));
    for (i=0; i < (NT-1); i++){
        T[i] = dt*((float)(i+1));
    }

    

    /*init SSAL */
    SSAL_Initialise(argc,argv);
    /*import model*/
    CRN = SSAL_ImportLSBML(filename);
    CRN_ptr = (SSAL_ChemicalReactionNetwork *)CRN.model;
    /*we will just use the default initial conditions*/
    X0 = CRN_ptr->X0;
    names = CRN_ptr->names;
    n = CRN_ptr->N;
    for(i=0;i<CRN_ptr->N;i++)
        fprintf(stderr,"%s\n",names[i]);
    /* build realisation simulation */
    switch (type)
    {
        default:
        case 1:
        {
            sim = SSAL_CreateRealisationsSim(&CRN,n,NULL,NR,NT,T,X0);
        }
            break;
        case 2:
        {
            sim = SSAL_CreateExpectedValueSim(&CRN,n,NULL,NR,NT,T,X0);
        }
            break;
    }

    free(T);
    SSAL_WriteChemicalReactionNetwork(stderr,*CRN_ptr);

    switch (algType)
    {
        default:
        case SSAL_ESSA_GILLESPIE_SEQUENTIAL:
            opts[0] = '\0';
            break;
        case SSAL_ASSA_TAU_LEAP_SEQUENTIAL:
            sprintf(opts,"--tau %f",tau);
            break;
        case SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL:
            sprintf(opts,"--tau0 %f --L %d --M %d --eps %f",tau0,L,M,eps);
            break;
    }
    /*simulate realisations*/
    SSAL_Simulate(&sim,algType,opts);
    /*write the data*/
    SSAL_WriteSimulation(stdout,sim);
    return 0;
}
