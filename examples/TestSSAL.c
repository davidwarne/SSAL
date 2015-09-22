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
#include<unistd.h>
#include "SSAL.h"

#define NT 5

int main(int argc , char ** argv)
{
    int i;
    int NR; /*number of realisations*/
    float dt;
    float T[NT];
    int model;
    char opts[256];
    char ch;
    SSAL_AlgorithmType algType;
    SSAL_Model CRN;
    SSAL_Simulation sim;

    NR = 1;
    model = 3;
    algType = SSAL_ESSA_GILLESPIE_SEQUENTIAL;
    T[NT-1] = 1.0;
    opts[0] = '\0';

    /*get input args with getops*/
    opterr = 0;
    while ((ch = getopt(argc,argv,"n:t:h:m:T:")) != -1)
    {
        switch(ch)
        {
            case 'n':
                NR = (int)atoi(optarg);
                break;
            case 't':
                T[NT-1] = (float)atof(optarg);
                break;
            case 'm':
                model = (int)atoi(optarg);
                break;
            case 'T': /*specifiying tau invokes tau leap method*/
                algType = SSAL_ASSA_TAU_LEAP_SEQUENTIAL;           
                sprintf(opts,"--tau %s",optarg);
                break;
            default:
            case '?':
                fprintf(stderr,"Unknown Option -%s\n",optarg);
            case 'h':
                fprintf(stderr,"Usage: %s [-n numRealisations] [-t endTime] [-m model] [-T tau]\n");
                break;
                break;
        }
    }
    
    /*build sample times*/
    dt = T[NT-1]/((float)(NT));
    for (i=0; i < (NT-1); i++){
        T[i] = dt*((float)(i+1));
    }

    /*init SSAL */
    SSAL_Initialise(argc,argv);
    /*select model*/
    switch(model)
    {
        default:
        case 1:
        {
            int m,n;
            m = 1;
            n = 1;
            {
                float nu_minus[1] = {1};
                float nu_plus[1] = {0};
                float c[1] = {0.1};
                char *names[1] = {"A"};
                float X0[1] = {20};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }
            break;
        case 2:
        {
            int m,n;
            m = 2;
            n = 1;
            {
                float nu_minus[2] = {1,0};
                float nu_plus[2] = {0,1};
                float c[2] = {0.1,1.0};
                char *names[1] = {"A"};
                float X0[1] = {0};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }

            break;
        case 3:
        {
            int n,m;
            m = 5;
            n = 3;
            {
                /*this is just the Dimerisation problem*/
                float nu_minus[15] = {0,0,0,1,0,0,0,2,0,1,0,0,0,1,0};
                float nu_plus[15] = {1,0,0,1,1,0,0,0,1,0,0,0,0,0,0};
                float c[5] = {25,1000,0.001,0.1,1};
                char *names[3] = {"M","P","D"};
                float X0[3] = {0,0,0};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }
            break;
        case 4:
        {
            int m,n;
            m = 4;
            n = 1;
            {
                float nu_minus[4] = {2,3,0,1};
                float nu_plus[4] = {3,2,1,0};
                float c[4] = {0.18,0.00025,2200,37.5};
                char *names[1] = {"A"};
                float X0[1] = {0};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }

            break;
        case 5:
        {
            int m,n;
            m = 4;
            n = 2;
            {
                float nu_minus[8] = {2,1,0,0,1,0,0,0};
                float nu_plus[8] = {3,0,1,0,0,0,0,1};
                float c[4] = {0.00004,50,10,25};
                char *names[2] ={"A","B"};
                float X0[2] = {10,10};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }
        

            break;
        case 6:
        {
            int m,n;
            m = 4;
            n = 3;
            {
                float nu_minus[12] = {1,0,0,0,1,0,2,0,0,0,1,0};
                float nu_plus[12] = {0,0,0,0,0,1,0,1,0,2,0,0};
                float c[4] = {0.04,0.002,0.5};
                char *names[3] = {"S1","S2","S3"};
                float X0[3] = {100000,0,0};
                /* build chemical reaction network*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,n,nu_minus,nu_plus,c);
                /* build realisation simulation */
                sim = SSAL_CreateRealisationsSim(
                    &CRN,n,names,NR,NT,T,X0);
            }
        }
       
            break;
    }

    /*write the CRN*/
    SSAL_WriteChemicalReactionNetwork(stdout,*((SSAL_ChemicalReactionNetwork *)(CRN.model)));
    /*simulate realisations*/
    SSAL_Simulate(&sim,algType,opts);
    /*write the data*/
    SSAL_WriteRealisationsSim(stdout,sim.sim);
    return 0;
}
