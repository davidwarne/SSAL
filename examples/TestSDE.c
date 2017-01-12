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
 * @file TestSDE.c
 * @brief An example program to generate correlated Euler-Maruyama realisations.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include<stdlib.h>
#include<stdio.h>
#include "SSAL.h"

/*example of SDEs using SSAL*/

/*OrnsteinUhlenbeck SDE*/

/* dX_t = -X_t dt + 0.3 dW_t*/

void 
mu(SSAL_real_t *X, unsigned int n, SSAL_real_t *params, unsigned int m, 
   SSAL_real_t t, SSAL_real_t* mur)
{
    mur[0] = -X[0];
}

void 
sigma(SSAL_real_t *X, unsigned int n, SSAL_real_t *params, unsigned int m, 
      SSAL_real_t t, SSAL_real_t* sigr)
{
    sigr[0] = params[0]; /*only one parameter*/
}



int main(int argc , char ** argv)
{
    int i,j;
    int NT;
    int NR; /*number of realisations*/
    SSAL_real_t dt;
    SSAL_real_t *T;
    SSAL_real_t *X0;
    SSAL_real_t *X_f_r;
    SSAL_real_t *X_c_r;
    SSAL_real_t sig;
    int d;
    char ** names;
    d = 0;
    //int model;
    int type;
    char opts[256];
    int m,n;
    sig = 0.3;
    names = (char **)malloc(sizeof(char *));
    names[0] = "X";

    X0 = (SSAL_real_t *)malloc(sizeof(SSAL_real_t));
    /*EM args*/
    SSAL_real_t h;
    SSAL_SDE SDE;
    NR = 3;
    NT = 1000;
    T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
    h = 5.0/((SSAL_real_t)NT);
    /*build sample times*/
    for (i=0; i < NT; i++){
        T[i] = h*((SSAL_real_t)(i+1));
    }

    
    X_f_r = (SSAL_real_t *)malloc(NR*NT*sizeof(SSAL_real_t));
    X_c_r = (SSAL_real_t *)malloc(NR*NT*sizeof(SSAL_real_t));

    /*init SSAL */
    SSAL_Initialise(argc,argv);
    /*build SDE model*/
    SDE = SSAL_CreateStochasticDifferentialEquation(names,1,1,&mu,&sigma,&sig);
    X0[0] = 2; 
    SDE.X0 = X0;
    names = SDE.names;
    n = SDE.N;
    for(i=0;i<SDE.N;i++)
        fprintf(stderr,"%s\n",names[i]);
    
    /*simulate realisations*/
    for (i=0;i<NR;i++)
    {
        dacems(1,1,NT,T,&sig,X0,&mu,&sigma,1,&d,h,2,X_f_r+i*NT,X_c_r+i*NT);    
        //daems(1,1,NT,T,&sig,X0,&mu,&sigma,1,&d,h,X_f_r+i*NT);    
    }

    /* Output sample paths*/
    fprintf(stdout,"\"R\",\"h\"");
    for (j=0;j<NT;j++)
    {
        fprintf(stdout,",\"T_%d\"",j);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"0,0");
    for (j=0;j<NT;j++)
    {
        fprintf(stdout,",%g",T[j]);
    }
    fprintf(stdout,"\n");
    for (i=0;i<NR;i++)
    {
        fprintf(stdout,"%d,%g",i,h);
        for (j=0;j<NT;j++)
        {
            fprintf(stdout,",%g",X_f_r[i*NT + j]);    
        }
        fprintf(stdout,"\n");
        fprintf(stdout,"%d,%g",i,h*2.0);
        for (j=0;j<NT;j++)
        {
            fprintf(stdout,",%g",X_c_r[i*NT + j]);    
        }
        fprintf(stdout,"\n");

    }
    return 0;
}
