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
 * @details approximates sampling of P(theta | D) by sampling theta then evaluating rho
 */
int ABCrejection(SSAL_Simulation *sim,int max_n, int nacc, float *data,
    int m,float *a, float *b, float (*rho)(int, float *, float *),
    float eps, float * theta, float *rhoVals, int *numAccept, 
    float *acceptRate)
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
        SSAL_Simulate(sim,SSAL_ESSA_GILLESPIE_SEQUENTIAL,NULL);

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
    float X0; /*initial condition*/
    float *a,*b; /*intervals of prior distributions*/
    float *c_real; /*rate parameters used to generate data*/
    float *X_data;
    char *names[1] = {"X"}; /*species symbol names*/
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
    X0 = 200.0;
    
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
            t_end = (float)atof(argv[++i]);
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
        else if (!strcmp("-S",argv[i]))
        {
            nt = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-X0",argv[i]))
        {
            X0 = (float)atof(argv[++i]);
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
            fprintf(stdout,"Usage: %s -N numSample -MaxN maxTrials -t t_end -m modelID M c1,...,cM -S timepoints -X0 initialCondition -e epsilon -a a1,...,aM -b b1,...bM\n",argv[0]);
            exit(1);
        }
    }

    

    /*init SSAL*/
    SSAL_Initialise(argc,argv);

    switch(model)
    {
        case 1:
        {
            int m;
            m = 1;
            {
                float nu_minus[1] = {1};
                float nu_plus[1] = {0};
                /*default parameter of 0.1 is a good start*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,1,nu_minus,nu_plus,c_real);
            }
        }            
            break;
        case 2:
        {
            int m;
            m = 2;
            {
                float nu_minus[2] = {1,0};
                float nu_plus[2] = {0,1};
                /*default parameters of {0.1, 1.0} is a good start*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,1,nu_minus,nu_plus,c_real);
            }
        } 
    }

    CRN_ptr = (SSAL_ChemicalReactionNetwork *)(CRN.model);


    /*allocate memory*/
    T = (float *)malloc(nt*sizeof(float));
    rhoV = (float *)malloc(nsamples*sizeof(float));
    c_sample = (float *)malloc(nsamples*CRN_ptr->M*sizeof(float));
     
    /*create our simulation*/
    T[nt -1] = t_end;
    for (i=0;i<(nt-1);i++)
    {
        T[i] = (T[nt-1])*(((float)i+1)/((float)nt));
    }
    sim = SSAL_CreateRealisationsSim(&CRN,1,names,1,nt,T,&X0);

    /*generate dummy data as an example*/
    //ExactSoln(model,nt,T, X0, CRN_ptr->c, X_data);

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
    
    /*output to file*/
    fprintf(stdout,"\"epsilon\",\"rho\",\"model\",\"AcceptRate\",\"numAccepts\"");
    for (j=0;j<M;j++)
    {
        fprintf(stdout,",\"c_real%d\"",j);
    }
    for (j=0;j<M;j++)
    {
        fprintf(stdout,",\"c_sample%d\"",j);
    }
    
    fprintf(stdout,"\n");
    if (numAccept != 0)
    {
        for (i=0;i<numAccept;i++)
        {
            fprintf(stdout,"%f,%f,%u,%f,%d",epsilon,sqrt(rhoV[i]),model,acceptRate,numAccept);
            for (j=0;j<M;j++)
            {
                fprintf(stdout,",%f",c_real[j]);
            }
            for (j=0;j<M;j++)
            {
                fprintf(stdout,",%f",c_sample[i*CRN_ptr->M+j]);
            }
            fprintf(stdout,"\n");
        }
    }
    else
    {
       fprintf(stdout,"%f,%f,%u,%f,%d",epsilon,0,model,0,numAccept);
       for (j=0;j<M;j++)
       {
           fprintf(stdout,",%f",c_real[j]);
       }
       for (j=0;j<M;j++)
       {
           fprintf(stdout,",%f",0);           
       }
       fprintf(stdout,"\n");
    }
}
