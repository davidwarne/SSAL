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
 * @file TestABC.c
 * @brief An example of using SSAL for Parameter inference using ABC-rejection
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensand University of Technology
 *
 * @date 1 Oct 2015
 */
#include<stdlib.h>
#include<stdio.h>
#include"SSAL.h"


/**
 * @brief distance measure
 * @details (||X - X*||/||X||)^2 where X* is the simulation data and 
 * X is the measured dataset.
 *
 * @param n dimensionality of data
 * @param X dataset 
 * @param X_star simulated data
 */
float rho(int n, float *X, float *X_star)
{
    int i;
    float d,sum; 
    sum = 0;

    for (i=0;i<n;i++)
    {
        sum += (X[i] - X_star[i])*(X[i]-X_star[i]);
    }

    d = sum;
    sum = 0;
    for (i=0;i<n;i++)
    {
        sum += X[i]*X[i];   
    }
    return (sum != 0) ? d/sum : d;
}

/**
 * @brief ABC-rejection scheme
 * @details approximates sampling of P(theta | D) by sampling \theta ~U(a,b) ,
 * then generating D* ~ s(\theta). The value of theta is accepted as a sample from
 *  the posterior if rho(D,D*) <= \eps
 *
 * @param sim SSAL_Simulation (currently only Reaslisations are supported)
 * @param alg the simulation algorithm to use
 * @param args algorithm argument string
 * @param max_n Maximum number of simulations to run
 * @param nacc the desired number of accepted simulations
 * @param data data set, assumed to be the average of measured data paths
 * @param m dimension of parameter space
 * @param a lower bounds of parameter distribution
 * @param b upper bounds of parameter distribution
 * @param rho The distance function for evaluating the similarity of D* to D
 * @param eps the threahold of sample acceptance
 * @param theta an array of samples accepted as samples of the posterior
 * @param rhoVals array that will contain the values of rho for accepted 
 * samples
 * @param numAccept the number of accepted samples collected
 * @param acceptRate numAccept/numTrials
 *
 */
int ABCrejection(SSAL_Simulation *sim, SSAL_AlgorithmType alg, char* args, int max_n, 
    int nacc, float *data, int m,float *a, float *b, float (*rho)(int, float *, float *),
    float eps, float * theta, float *rhoVals, int *numAccept, float *acceptRate)
{
    int i,j,k;
    SSAL_ChemicalReactionNetwork *CRN_ptr;
    SSAL_RealisationSimulation *RS_ptr;
    float *X_r;

    CRN_ptr = (SSAL_ChemicalReactionNetwork *)(sim->model->model);
    RS_ptr = (SSAL_RealisationSimulation *)(sim->sim);
    X_r = (float *)(RS_ptr->output);
    
    k=0;
    for (i=0;i<max_n;i++)
    {
        float d;
        /*generate sample theta ~ U(a,b)*/
        for (j=0;j<m;j++)
        {
            CRN_ptr->c[j] = surngus(a[j],b[j]);     
        }
        
        /*run simulation*/
        SSAL_Simulate(sim,alg,args);

        /*accept/reject*/
        d = (*rho)((RS_ptr->Nvar)*(RS_ptr->NT),data,X_r);
        if (d <= eps*eps)
        {
            rhoVals[k] = d;
            for (j=0;j<m;j++)
            {
                theta[k*m + j] = CRN_ptr->c[j]; 
            }
            k++;

            if (k == nacc)
            {
                break;
            }
        }
    }

    *numAccept = k;
    *acceptRate = ((float)k)/((float)i);

    return 0;
}


/**
 * @brief Exact solution for deterministic versions of model
 * 
 */
int ExactSoln(int model,int nt, float * T, float X0, float *c, float *X)
{
    int t;
    if (model == 1)
    {
        for (t=0;t<nt;t++)
        {
            X[t] = X0*expf(-c[0]*T[t]);
        }
    }
    else
    {
        for (t=0;t<nt;t++)
        {
            X[t] = c[1]/c[0] + (X0 - c[1]/c[0])*expf(-c[0]*T[t]);
        }
    }
    return 0;
}

int WriteData(FILE *stream, int M,int numAccept,float epsilon, float *rhoV,int model,
    float acceptRate,float * theta_r, float *theta_s)
{
    static unsigned char header = 1;

    /*output to file*/
    if (header == 1)
    {
        fprintf(stream,"\"model\",\"epsilon\",\"AcceptRate\",\"numAccepts\"");
        for (j=0;j<M;j++)
        {
            fprintf(stream,",\"theta_r%d\"",j);
        }
        fprintf(",\"rho\"");
        for (j=0;j<M;j++)
        {
            fprintf(stream,",\"theta_s%d\"",j);
        }
        fprintf(stream,"\n");
        header = 0;
    }
 
    for (i=0;i<numAccept;i++)
    {
        fprintf(stream,"%u,%f,%f,%d",model,epsilon,acceptRate,numAccept);
        for (j=0;j<M;j++)
        {
            fprintf(stream,",%f",theta_r[j]);
        }
        fprintf(stream,",%f",rhoV[i]);
        for (j=0;j<M;j++)
        {
            fprintf(stream,",%f",theta_s[i*M+j]);
        }
        fprintf(stream,"\n");
    }
 }

int main(int argc,char ** argv)
{
    float epsilon; /*acceptance threshold*/
    float *rhoV; /*array to store distances*/
    float *c_sample; /*array to store samples*/
    int numAccept;
    float acceptRate;
    unsigned int genData;

    int nsamples; /*number of samples from posterior to obtain*/
    int nMax; /*Max number of simulations runs*/
    float *T; /*time points as out summary statistic*/
    float t_end; /*simulation endtim*/
    int nt; /*size of summary statistic*/
    int model; /*1 = degradation 2 = production/degradation*/
    int M; /*number of parameters specified*/
    float *X0; /*initial condition*/
    int N; /*size of state space*/
    float *a,*b; /*intervals of prior distributions*/
    float *c_real; /*rate parameters used to generate data*/
    float *X_data;
    char **names; /*species symbol names*/
    int i,j; /*loop counters*/

    SSAL_Model CRN; /*Chemical reaction model*/
    SSAL_ChemicalReactionNetwork *CRN_ptr;
    SSAL_ExpectedValueSimulation *EV_ptr;
    SSAL_Simulation sim; /*simulation*/
    SSAL_Simulation simData; /*simulation*/
    
    /*default values */
    nt = 1;
    t_end = 30.0;
    model = 1;
    M = 1;
    
    nsamples = 10000;
    nMax = 1000000;
    genData = 0;
    /*get input args*/
    for (i=1;i<argc;i++)
    {
        
        if(!strcmp("-N",argv[i]))
        {
            nsamples = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-d",argv[i]))
        {
            genData = 1;
        }
        else if (!strcmp("-MaxN",argv[i]))
        {
            nMax = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-t",argv[i]))
        {
            nt = (int)atoi(argv[++i]);
            T = (float *)malloc(nt*sizeof(float));
            for (j=0;j<NT;j++)
            {
                T[j] = (float)atof(argv[++i]);
            }
        }
        else if (!strcmp("-m",argv[i]))
        {
            model = (int)atoi(argv[++i]);
            M = (int)atoi(argv[++i]);
            c_real = (float *)malloc(M*sizeof(float));
            for (j=0;j<M;j++)
            {
                c_real[j] = (float)atof(argv[++i]);
            }
        }
        else if (!strcmp("-X0",argv[i]))
        {
            N = (int)atoi(argv[++i]);
            X0 = (float *)malloc(N*sizeof(float));
            for (j=0;j<N;j++)
            {
                X0[j] = (float)atof(argv[++i]);
            }
        }
        else if (!strcmp("-e",argv[i]))
        {
            epsilon = (float)atof(argv[++i]);
        }
        else if (!strcmp("-a",argv[i]))
        {
            a = (float *)malloc(M*sizeof(float));
            for (j=0;j<M;j++)
            {
                a[j] = (float)atof(argv[++i]);
            }
        }
        else if (!strcmp("-b",argv[i]))
        {
            b = (float *)malloc(M*sizeof(float));
            for (j=0;j<M;j++)
            {
                b[j] = (float)atof(argv[++i]);
            }
        }
        else
        {
            fprintf(stdout,"Usage: %s -N numAccepts -MaxN maxTrials -m modelID M c1,...,cM -t NT t1,..,tNT -X0 N X_1(0),...,X_N(0)  -e epsilon -a a1,...,aM -b b1,...bM -T tau\n",argv[0]);
            exit(1);
        }
    }

    /*allocate memory of reaction network*/
    names = (char **)malloc(N*sizeof(char *));
    nu_minus = (float *)malloc(N*M*sizeof(float));
    nu_plus = (float *)malloc(N*M*sizeof(float));

    for (i=0;i<N*N:i++)
    {
        nu_minus[i] = 0;; 
    }
    for (i=0;i<N*N:i++)
    {
        nu_plus[i] = 0;; 
    }

    /*init SSAL*/
    SSAL_Initialise(argc,argv);

    switch(model)
    {
        case 1: /*simple degradation*/
            names[0] = "A";
            nu_minus[0] = 1;
            break;
        case 2: /*production/degradation*/ 
            /*note: this is actually redundant with the monomolecular chain*/
            names[0] = "A";
            nu_minus[0] = 1; /*A --> 0*/
            nu_plus[1] = 1; /*0 --> A*/
            break;
        case 3: /*monomolecular chain 0 --> X_1 --> ... X_i-1 --> X_i --> ... X_N --> 0*/
            for (i=0;i<N;i++)
            {
                names[i] = (char *)malloc(SSAL_MAX_NAME_SIZE*sizeof(char));
                sprintf(names[i],"X_%d",i);
            }

            
            nu_plus[0] = 1; /*0 --> X_1*/
            for (i=1;i<(N-1);i++)
            {
                nu_minus[i*N + (i -1)] = 1;       
                nu_plus[i*N + i] = 1;   /*X_{i-1} --> X_{i}, i = 2,..,N-1*/ 
            }
            nu_minus[N*N - 1] = 1; /*X_N --> 0*/
            break;
        case 4: /*random monomolecular reaction network*/
            /* conversion: X_i --> X_j i != j
             * production: 0 --> X_i
             * degradation: X_i --> 0
             */
            for (i=0;i<N;i++)
            {
                names[i] = (char *)malloc(SSAL_MAX_NAME_SIZE*sizeof(char));
                sprintf(names[i],"X_%d",i);
            }
            
            for (j=0;j<M;j++)
            {
                if (rand() > (RAND_MAX >> 2)) /*add conversion reaction*/
                {
                    int ii,jj;
                    ii = rand();
                    jj = rand();

                    nu_minus[j*N + ii%N] = 1;
                    if (ii != jj)
                    {
                        nu_plus[j*N + jj%N] = 1;
                    }
                    else
                    {
                        nu_plus[j*N + (jj+1)%N];
                    }
                }
                else
                {
                    if (rand() > (RAND_MAX >> 1)) /*add production reaction*/
                    {
                        nu_plus[j*N + rand()%N] = 1;
                    }
                    else /*add degradation reaction*/
                    {
                        nu_minus[j*N + rand()%N] = 1;
                    }
                }
            }
            break;
    }

    /*create chemical reaction network*/
    CRN = SSAL_CreateChemicalReactionNetwork(names,M,N,nu_minus,nu_plus,c_real);
    CRN_ptr = (SSAL_ChemicalReactionNetwork *)(CRN.model);


    /*allocate memory*/
    rhoV = (float *)malloc(nsamples*sizeof(float));
    c_sample = (float *)malloc(nsamples*CRN_ptr->M*sizeof(float));
     
    
    sim = SSAL_CreateRealisationsSim(&CRN,1,names,1,nt,T,&X0);

    /*generate dummy data as an example*/
    simData = SSAL_CreateExpectedValueSim(&CRN,1,names,100,nt,T,&X0);
    EV_ptr = (SSAL_ExpectedValueSimulation *)(simData.sim); 
    SSAL_Simulate(&simData,SSAL_ESSA_GILLESPIE_SEQUENTIAL,NULL);
    if (genData)
    {
        SSAL_WriteSimulation(stdout,simData);
        exit(0);
    }

    X_data = EV_ptr->E;
    for (i=0;i<nt;i++)
    {
        X_data[i] = round(X_data[i]);
    }

    /*apply ABC rejection method*/
    ABCrejection(&sim,nMax,nsamples, X_data,CRN_ptr->M, a, b, rho, 
            epsilon,c_sample, rhoV, &numAccept, &acceptRate);
    
    /*write results*/
    WriteData(stdout, CRN_ptr->M, numAccept, epsilon, rhoV, model,
        acceptRate,c_real, c_sample);
    exit(0);
}
