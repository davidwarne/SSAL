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
float rho(n,X,X_star)
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
    return d/sum;
}

/**
 * @brief ABC-rejection scheme
 * @details approximates sampling of P(theta | D) via 
 */
int ABCrejection(SSAL_Simulation *sim,int n, float *data,int m,float *a, float *b, float * theta, float *rho)
{
    int i;
    for (i=0;i<nsamples;i++)
    {
        /*generate sample in theta*/
        for (j=0;j<m;j++)
        {
            
        }
        
        SSAL_Simulate
    }
    return 0;
}


/**
 * @brief Exact solution for deterministic versions of model
 * 
 */
int ExactSoln(int model,int nt, float * T, float X0, float *c, float *X)
{
    int t;
    if (m==1)
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
            X[t] = c[0]/c[1] + (X0 - c[0]/c[1])*expf(-c[1]*T[t]);
        }
    }
    return 0;
}

int main(int argc,char ** argv)
{
    float epsilon; /*acceptance threshold*/
    float *rho; /*array to store distances*/
    float *c_sample; /*array to store samples*/

    int nsamples; /*number of samples to trial*/
    float *T; /*time points as out summary statistic*/
    float t_end; /*simulation endtim*/
    int nt; /*size of summary statistic*/
    int model; /*1 = degradation 2 = production/degradation*/
    int M; /*number of parameters specified*/
    float X0; /*initial condition*/
    float *a,*b; /*intervals of prior distributions*/
    float *c_real; /*rate parameters used to generate data*/
    char *names[1] = {"X"}; /*species symbol names*/
    int i,j; /*loop counters*/
    SSAL_Model CRN; /*Chemical reaction model*/
    SSAL_Simulation sim; /*simulation*/
    
    /*default values */
    nt = 1;
    t_end = 30.0;
    model = 1
    M = 1;
    X0 = 200.0;
    
    nsamples = 10000;
    /*get input args*/
    for (i=0;i<argc;i++)
    {
        if(!strcmp("-N",argv[i]))
        {
            nsamples = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-t",argv[i]))
        {
            t_end = (float)argv[++i];
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
            epsilon = (float)argv[++i];
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
                /*default parameters of {0.1, 1.0} is a sgood start*/
                CRN = SSAL_CreateChemicalReactionNetwork(
                    names,m,1,nu_minus,nu_plus,c_real);
            }
        } 
    }


    /*allocate memory*/
    X_data = (float *)malloc(nt*sizeof(float));
    T = (float *)malloc(nt*sizeof(float));
    rho = (float *)malloc(nsamples*sizeof(float));
    c_sample = (float *)malloc(nsamples*CRN.M*sizeof(float));
 
    /*create our simulation*/
    T[nt -1] = t_end;
    for (i=0;i<(nt-1);i++)
    {
        T[i] = (T[nt-1])*(((float)i+1)/((float)nt));
    }
    sim = SSAL_CreateRealisationsSim(&CRN,1,names,1,nt,T,X0);
    
    /*generate dummy data as an example*/
    ExactSoln(model,nt,T, X0, CRN.c, X_data);

    /*apply ABC rejection method*/
    ABCrejection(&sim,,nsamples,X_data,a,b,c_sample,rho);

    /*output to file*/
    fprintf(stdout,"\"epsilon\",\"rho\",\"accepted\"");
    for (j=0;j<M;j++)
    {
        fprintf(stdout,",\"c_real%d\"",j);
    }
    for (j=0;j<M;j++)
    {
        fprintf(stdout,",\"c_sample%d\"",j);
    }
    fprintf(stdout,"\n");
    for (i=0;i<nsamples;i++)
    {
        fprintf(stdout,"%f,%f,%u",epsilon,sqrt(rho[i]),(rho[i] < epsilon*epsilon));
        for (j=0;j<M;j++)
        {
            fprintf(stdout,",%f",c_real[j]);
        }
        for (j=0;j<M;j++)
        {
            fprintf(stdout,",%f",c_sample[i*CRN.M+j]);
        }
        fprintf(stdout,"\n");
    }
}
