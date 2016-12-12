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

/*example of ODEs using RK4*/

/*Van der Pol ODE*/

/* dx/dt = y , dy/dt = \mu(1-x^2)y - x*/

void vdp(SSAL_real_t *Y, unsigned int n, SSAL_real_t *params, unsigned int m, SSAL_real_t t, SSAL_real_t* f_r){
    f_r[0] = Y[1];
    f_r[1] = params[0]*(1.0 - Y[0]*Y[0])*Y[1] - Y[0];
}

/*Lotka-Volterra*/
/* dx/dt = ax - yx , dy/dt = byx - y*/
void lv(SSAL_real_t *Y, unsigned int n, SSAL_real_t *params, unsigned int m, SSAL_real_t t, SSAL_real_t* f_r){
    f_r[0] = params[0]*Y[0] - Y[0]*Y[1];
    f_r[1] = params[1]*Y[0]*Y[1] - Y[1];
}



int main(int argc , char ** argv)
{
    int i;
    int NT;
    SSAL_real_t dt;
    SSAL_real_t *T;
    SSAL_real_t *Y0;
    SSAL_real_t *Y_r;
    SSAL_real_t *params;
    SSAL_real_t h;
    int *d;
    int type;
    char opts[256];
    int m,n;

    SSAL_Initialise(argc,argv);
    //m = 1;
    m = 2;
    n = 2;
    NT = (int)atoi(argv[1]);

    Y0 = (SSAL_real_t *)malloc(n*sizeof(SSAL_real_t));
    params= (SSAL_real_t *)malloc(m*sizeof(SSAL_real_t));
    d = (int *)malloc(n*sizeof(int));
    T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
    Y_r = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
    
    //h = 100.0/((SSAL_real_t)NT);
    h = 15.0/((SSAL_real_t)NT);
    /*build sample times*/
    for (i=0; i < NT; i++){
        T[i] = h*((SSAL_real_t)(i+1));
    }
    // use small timestep internally
    h = 15.0/10000.0;
    d[0] = 0;
    d[1] = 1;
    params[0]= 1.0;
//    params[0]= (double)atof(argv[1]);
    params[1]= 1.0;
    Y0[0] = 1.0; 
    Y0[1] = 0.5; 
    //Y0[0] = 1e-6; 
    //Y0[1] = 1e-6; 
    drk4s(m,n,NT,T,params,Y0,&lv,n,d,h,Y_r);

    for (i=0;i<NT;i++)
    {
        // peturb with noise
    //    fprintf(stdout,"%f,%f,%f\n", T[i],Y_r[i]+durngns(0.0,0.5),Y_r[NT+i]+ durngns(0.0,0.5));
        fprintf(stdout,"%f,%f,%f\n", T[i],Y_r[i],Y_r[NT+i]);
    }
    return 0;
}
