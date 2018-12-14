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

#ifndef __SSAL_H_
#define __SSAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "cJSON.h"
#include "SSAL_internal.h"


/*error codes*/
#define SSAL_SUCCESS                    0
#define SSAL_INVALID_OPTION_ERROR       -1
#define SSAL_UNKNOWN_OPTION_ERROR       -2
#define SSAL_IO_ERROR                   -3
#define SSAL_UNKNOWN_TYPE_ERROR         -4
#define SSAL_MEMORY_ERROR               -5
#define SSAL_UNSUPPORTED_ERROR          -6

/*option codes*/
#define SSAL_MAX_NAME_SIZE      128
#define SSAL_MAX_BUFFER_SIZE    1024


typedef double SSAL_real_t;

/**
 * @struct SSAL_ChemicalReactionNetwork_struct
 * @brief Defines a stochastic Chemical Reaction network model
 * @detail A user can manually interact with an instance of this structure
 * or use the SSAL_Create_ChemicalReactionNetwork function
 */
struct SSAL_ChemicalReactionNetwork_struct 
{
    /** Chemical Species Names*/
    char ** names;
    /** number of Chemical Species*/
    uint32_t N;
    /** number of reactions */
    uint32_t M;
    /** stochiometric coefficient matrix of the reactants */
    SSAL_real_t *nu_minus;
    /** stochiometric coefficient matrix of the products */
    SSAL_real_t *nu_plus;
    /** stochiometric matrix nu =  (nu^+ -nu^-) */
    SSAL_real_t *nu;
    /** kinetic rate constants*/
    SSAL_real_t *c;
    /**default initial condition*/
    SSAL_real_t *X0;
    /** number of variables to track */
    int nvar;
    /** designates which variables to track */
    int *vars;
};

/** type name for  SSAL_ChemicalReactionNetwork_struct*/
typedef struct SSAL_ChemicalReactionNetwork_struct SSAL_CRN;

/**
 * @struct SSAL_StochasticDifferentialEquation_struct
 * @brief general vector SDE of the form
 *          dX_t = mu(X_t,t)dt + sigma(X_t,t)dW_t 
 * where X_t \in R^N and W_t is a vector of N Wiener processes.
 */
struct SSAL_StochasticDifferentialEquation_struct 
{
    /** Variable Names*/
    char ** names;
    /** number of equations */
    uint32_t N;
    /** number of parameters */
    uint32_t M;
    /** drift function*/
    void (*mu)(SSAL_real_t*,uint32_t, SSAL_real_t*, uint32_t, 
               SSAL_real_t, SSAL_real_t*);
    /** diffusion function */
    void (*sigma)(SSAL_real_t*,uint32_t,SSAL_real_t *, uint32_t, 
                  SSAL_real_t,SSAL_real_t*);
    /**@todo add corvariance matrix support for true multi-variate SDE
     * currently only SDE's driven by independent standard Brownian motions 
     * are supported
     */    
    /** parameter vector */
    SSAL_real_t *p;
    /** default initial conditions*/
    SSAL_real_t *X0;
};
typedef struct SSAL_StochasticDifferentialEquation_struct SSAL_SDE;


/*macro functions*/


/*converts variable name list to indexes*/
#define SSAL_VARS2INDS(Nvar,var,model,varInd) {             \
    int i,j;                                                \
    if (Nvar == (model)->N)                                 \
    {                                                       \
        for (j=0;j<(model)->N;j++)                          \
        {                                                   \
            varInd[j] = j;                                  \
        }                                                   \
        for (i=0;i<Nvar;i++)                                \
        {                                                         \
            strncpy(var[i],(model)->names[i],SSAL_MAX_NAME_SIZE); \
        }                                                         \
    }                                                       \
    for (i=0;i<Nvar;i++)                                    \
    {                                                       \
        for (j=0;j<(model)->N;j++)                          \
        {                                                   \
            if (!strcmp(var[i],(model)->names[j]))          \
            {                                               \
                varInd[i] = j;                              \
            }                                               \
        }                                                   \
    }                                                       \
}

/* function prototypes*/

/*==================================================================*/
/*           SSAL Generic API Functions                             */
/*==================================================================*/

/*initialisation an error handling*/
void 
SSAL_InitRNGS(unsigned int *, int);

void 
SSAL_HandleError(int , char *,int, unsigned char, unsigned char, char *); 

int 
SSAL_Initialise(int,char **); 

/* generic wrapper functions*/
SSAL_CRN 
SSAL_ImportLSBML(const char * );

/*==================================================================*/
/*         Chemical Reaction Network API                            */
/*==================================================================*/

/*object creation*/
SSAL_CRN 
SSAL_CreateChemicalReactionNetwork(char **, int ,int, SSAL_real_t * restrict , 
                                   SSAL_real_t * restrict , 
                                   SSAL_real_t * restrict );

int 
SSAL_WriteChemicalReactionNetwork(FILE *,SSAL_CRN);

/*==================================================================*/
/*         Stochastic Differential Equation API                     */
/*==================================================================*/
/* SDE structures and functions*/
SSAL_SDE 
SSAL_CreateStochasticDifferentialEquation(char ** , int, int ,
                                          void (*)(SSAL_real_t *,uint32_t,
                                                   SSAL_real_t*, uint32_t, 
                                                   SSAL_real_t,SSAL_real_t*),
                                          void (*)(SSAL_real_t *,uint32_t,
                                                   SSAL_real_t*,uint32_t, 
                                                   SSAL_real_t, SSAL_real_t*), 
                                          SSAL_real_t *);

/*utilities*/
char** 
SSAL_UtilTokeniseArgs(int *,const char *);
#endif /*__SSAL_H_ */
