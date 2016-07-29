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
    int d; /* dimesion of random vector*/
    double mu[2];
    double lambda[4];
    double Z[2];
    double X[2];

    NR = 10000;
    d = 2;
    mu[0] = 1.0;
    mu[1] = 0.5;
    lambda[0] = 1.0;
    lambda[1] = 0.0;
    lambda[2] = 0.3;
    lambda[3] = 0.7;
    

    /*init SSAL */
    SSAL_Initialise(argc,argv);

    for (i=0;i<NR;i++)
    {
        durngmvns(d,mu,lambda,Z,X);
        fprintf(stdout,"%f,%f,%f,%f\n",Z[0],Z[1],X[0],X[1]);
    }
    return 0;
}
