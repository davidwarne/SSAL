# SSAL
A fast Stochastic Simulation Algorithm Library

## Summary

A number of routines for numerical simulation of stochastic processes. This project is very much research code, and is still under development.

## Developer
David J. Warne (david.warne@qut.edu.au), School of Mathematical Sciences, Science and Engineering Faculty, Queensland University of Technology.


## License

SSAL: Stochastic Simulation Algorithm Library
Copyright (C) 2018  David J. Warne

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Current Features

1. Varitety of routines for sampling standard probability distributions.
2. Exact Stochastic Simulation Algorithms (ESSA) for discrete state, continuous-time Markov processes
3. Approximate Stochastic Simulatio Algorithms (ASSA) for discrete state, continuous-time Markov processes and Stochastic Differential Equations (SDEs), including coupled schemes ideal for multilevel Monte Carlo (MLMC). 

Currently, the library is all sequential with minimal library dependences (the Intel MKL_VSL library is optional to imporve RNG performance).

The library has been designed with performance and usability in mind. 

## Future development plans

Many of these rountines I use in my own research. Therefore, the current features and level of performance etc will generally be updated based on my own usage of the library. I am planning parallel (multi-threaded and distributed) version of many of the rountines.

If you find this library useful, or have any suggestions for improvement, I am always happy to be contacted.

## Example usage

A number of examples are provide in the `examples/` folder. However, here is a minimalistic example for SDE simulation using the Euler-Maruyama scheme

```C
#include<stdlib.h>
#include<stdio.h>
#include "SSAL.h"

/*example of SDEs using SSAL*/
/*OrnsteinUhlenbeck SDE dX_t = -X_t dt + 0.3 dW_t*/
/* drift function */
void 
mu(SSAL_real_t *X, unsigned int n, SSAL_real_t *params, unsigned int m, 
   SSAL_real_t t, SSAL_real_t* mur) 
{
    mur[0] = -X[0];
}

/* diffusion function*/
void 
sigma(SSAL_real_t *X, unsigned int n, SSAL_real_t *params, unsigned int m, 
      SSAL_real_t t, SSAL_real_t* sigr)
{
    sigr[0] = params[0]; /*only one parameter*/
}

int 
main(int argc , char ** argv)
{
    int NT, NR;                  /* number of observation times and realisations*/
    SSAL_real_t sig;             /* diffusion coefficient*/
    SSAL_real_t h;               /* time step*/
    SSAL_real_t *T, *X0, *X_f_r; /* store times, initial condition and realisations*/
    char **names;                /* array of SDE variable names*/
    int d[1];                    /* array of SDE variable indices*/
    SSAL_SDE SDE;                /* SDE structure*/
    
    /* initialise variables*/
    NT = 1000; NR = 3;
    sig = 0.3; h = 5.0/((SSAL_real_t)NT);  
    T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));
    /*build sample times*/
    for (int i=0; i < NT; i++) {
        T[i] = h*((SSAL_real_t)(i+1));
    }
    X0 = (SSAL_real_t *)malloc(sizeof(SSAL_real_t));
    X_f_r = (SSAL_real_t *)malloc(NR*NT*sizeof(SSAL_real_t)); 
    names = (char **)malloc(sizeof(char *));
    names[0] = "X";
    d[0] = 0;
  
    /*init SSAL */
    SSAL_Initialise(argc,argv);
    
    /*build SDE model*/
    SDE = SSAL_CreateStochasticDifferentialEquation(names,1,1,&mu,&sigma,&sig);
    X0[0] = 2; SDE.X0 = X0; names = SDE.names; n = SDE.N;
    
    /*simulate realisations*/
    for (int i=0;i<NR;i++) {
        daems(1,1,NT,T,&sig,X0,&mu,&sigma,1,&d,h,X_f_r+i*NT);    
    }

    /* Output final time steps*/
    fprintf(stdout,"\"Realisation\",\"X[T]\"");
    for (int i=0;i<NR;i++) {
        fprintf(stdout,"%d,%g\n",i,X_f_r[(i+1)*NT-1]);
    }
    return 0;
}

```

