/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2018  David J. Warne
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
 * @file TestPhaseTypeDist.c
 * @brief Example program for estimating the Phase-type distribution of a 
 * stochastic reaction-diffusion model.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include "SSAL.h"

/**
 * @brief Creates 1d instance of reaction-diffusion model. 
 * @details Let X be a chemical species and X_i(t) be the copy number of X at 
 * time t in the ith volume of dimension h, i.e., the copy number in 
 * [ih,(i+1)h). Let L = Nh and i = 0,1,...,N-1. 
 *
 *           k         k         k         k         k
 *     X_0 <---> X_1 <---> X_2 <---> ... <---> X_N-1 --> A
 * |<---h--->|<---h--->|<---h--->|   ...   |<---h--->)
 * where k = D/h^2 where D is the diffusivity in the continuum limit and A is 
 * an absorptive boundary. X_0 implicitly includes a reflective boundary
 *
 * @todo would like to add a Reaction Diffusion structure to simplify this 
 * long term a bit. could have: 1) a CRN; 2) diffusivity function;
 * and 3) boundaries/interfaces. The some functions to convert to compartment-based
 * RDME or SDE/indiviual based model.
 */
SSAL_CRN CreateSpatialCRN(SSAL_real_t L, SSAL_real_t h, SSAL_real_t D)
{
    int N; /*number of cells*/
    SSAL_CRN newCRN;
    char ** names;          /*strings of printable variable names*/
    SSAL_real_t * nu_minus; /* nu^- reactant stoichiometric coefficients*/
    SSAL_real_t * nu_plus;  /* nu^+ product stoichiomeric coefficients*/
    SSAL_real_t * c;        /* rate parameters c[i] = k for all i here*/

    N = (int)ceil(L/h); /* adjust h if needed to ensure L = Nh and N is int*/
    if (L < ((SSAL_real_t)N)*h)
    {
        h = L/((SSAL_real_t)N);
    }

    /* build arrays required for the CRN*/
    names = (char**)malloc((N+1)*sizeof(char*));
    names[0] = (char*)malloc((N+1)*SSAL_MAX_NAME_SIZE*sizeof(char));
    for (int i=1;i<N+1;i++)
    {
        names[i] = names[i-1] + SSAL_MAX_NAME_SIZE;
    }
    nu_minus = (SSAL_real_t *) malloc((N+1)*(2*N -1)*sizeof(SSAL_real_t));
    nu_plus = (SSAL_real_t *) malloc((N+1)*(2*N -1)*sizeof(SSAL_real_t));
    c = (SSAL_real_t *) malloc((2*N -1)*sizeof(SSAL_real_t));
    
    /* set names */
    for (int i=0;i<N;i++)
    {
        snprintf(names[i],SSAL_MAX_NAME_SIZE,"X_%d",i+1); 
    }
    snprintf(names[N],SSAL_MAX_NAME_SIZE,"A");
    
    /* set stoichiometries*/
    memset(nu_minus,0,(N+1)*(2*N+1)*sizeof(SSAL_real_t));
    memset(nu_plus,0,(N+1)*(2*N+1)*sizeof(SSAL_real_t));
    nu_minus[0] = 1.0;
    nu_plus[1] = 1.0;
    for (int i=1;i<N;i++)
    {
        /* set reactant stoichiometries*/
        nu_minus[(2*i-1)*(N+1) + i] = 1.0; 
        nu_minus[(2*i)*(N+1) + i] = 1.0;
        /* set product stoichiometries*/
        nu_plus[(2*i-1)*(N+1) + i -1] = 1.0; 
        nu_plus[(2*i)*(N+1) + i + 1] = 1.0;
    }
    /*rate parameters*/
    for (int i=0;i<2*N-1;i++)
    {
        /* compute rate constant to ensure diffusivity D is recovered for h -> 0*/
        c[i] = D/(h*h);
    }

    /* construct the CRN*/
    newCRN = CreateChemicalReactionNetwork(names,2*N-1,N+1,nu_minus,nu_plus,c);
    free(names[0]);
    free(names);
    free(nu_minus);
    free(nu_plus);
    free(c);
}


int
main(int argc, char ** argv)
