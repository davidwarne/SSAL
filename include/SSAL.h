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
#define SSAL_UNSUPPORTED_ALGORITHM_ERROR -7

/*option codes*/
#define SSAL_MAX_NAME_SIZE      128
#define SSAL_MAX_BUFFER_SIZE    1024


typedef double SSAL_real_t;

/**
 * @enum SSAL_AlgorithmType_enum
 * @brief Enumeration of the available Exact and Approximate algorithms
 */
enum SSAL_AlgorithmType_enum {
    SSAL_ESSA_GILLESPIE_SEQUENTIAL,
    SSAL_ESSA_GIBSON_BRUCK_SEQUENTIAL,
    SSAL_ASSA_TAU_LEAP_SEQUENTIAL,
    SSAL_ASSA_TAU_LEAP_CORRELATED_SEQUENTIAL,
    SSAL_ASSA_EULER_MARUYAMA_SEQUENTIAL,
    SSAL_ASSA_EULER_MARUYAMA_CORRELATED_SEQUENTIAL
};
/** type name for  SSAL_Algorithm_enum*/
typedef enum SSAL_AlgorithmType_enum SSAL_AlgorithmType;

enum SSAL_ModelType_enum {
    SSAL_CHEMICAL_REACTION_NETWORK,
    SSAL_STOCHASTIC_DIFFERENTIAL_EQUATION
};
/** type name for  SSAL_ModelType_enum*/
typedef enum SSAL_ModelType_enum SSAL_ModelType;

enum SSAL_SimulationType_enum {
    SSAL_REALISATIONS
};
/** type name for  SSAL_SimulationType_enum*/
typedef enum SSAL_SimulationType_enum SSAL_SimulationType;

/**
 * generic stochastic model 
 */
struct SSAL_Model_struct {
    SSAL_ModelType type;
    void *model;
};
/** type name for  SSAL_Model struct*/
typedef struct SSAL_Model_struct SSAL_Model;

/**
 * @struct SSAL_ChemicalReactionNetwork_struct
 * @brief Defines a stochastic Chemical Reaction network model
 * @detail A user can manually interact with an instance of this structure
 * or use the SSAL_Create_ChemicalReactionNetwork function
 */
struct SSAL_ChemicalReactionNetwork_struct {
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
    /** kinetic rate constants*/
    SSAL_real_t *c;
    /**default initial condition*/
    SSAL_real_t *X0;
};
/** type name for  SSAL_ChemicalReactionNetwork_struct*/
typedef struct SSAL_ChemicalReactionNetwork_struct SSAL_ChemicalReactionNetwork;

/**
 * @struct SSAL_StochasticDifferentialEquation_struct
 * @brief general vector SDE of the form
 *          dX_t = mu(X_t,t)dt + sigma(X_t,t)dW_t 
 * where X_t \in R^N and W_t is a vector of N Wiener processes.
 */
struct SSAL_StochasticDifferentialEquation_struct {
    /** Variable Names*/
    char ** names;
    /** number of equations */
    uint32_t N;
    /** number of parameters */
    uint32_t M;
    /** drift function*/
    void (*mu)(SSAL_real_t*,uint32_t, SSAL_real_t*, uint32_t, SSAL_real_t, SSAL_real_t*);
    /** diffusion function */
    void (*sigma)(SSAL_real_t*,uint32_t,SSAL_real_t *, uint32_t, SSAL_real_t,SSAL_real_t*);
    /**@todo add corvariance matrix support for true multi-variate SDE
     * currently only SDE's driven by independent standard Brownian motions are supported
     */    
    /** parameter vector */
    SSAL_real_t *p;
    /** default initial conditions*/
    SSAL_real_t *X0;
};
typedef struct SSAL_StochasticDifferentialEquation_struct SSAL_StochasticDifferentialEquation;

/**
 * @struct SSAL_Simulation_struct
 * @brief defines generic simulation 
 */
struct SSAL_Simulation_struct{
    /**
     * The type of simulation structure
     */
    SSAL_SimulationType type;
    /** model to use in simulation */
    SSAL_Model *model;
    void *sim;
};
/** type name for  SSAL_Simulation_struct*/
typedef struct SSAL_Simulation_struct SSAL_Simulation;

/**
 * realisation simulation
 */
struct SSAL_RealisationSimulation_struct {
    /**numbr of observation times */
    int NT;
    /** observation times */
    SSAL_real_t *T;
    /** number of realisations*/
    int NR;
    /**number of initial conditions */
    int NIC;
    /** initial conditions */
    SSAL_real_t *IC;
    /** number of variables to observe */
    int Nvar;
    /** model variables to observe by name */
    char **var;
    int *varInd;
    /** output data*/
    SSAL_real_t *output;
};
/** type name for  SSAL_RealisationSimulation_struct*/
typedef struct SSAL_RealisationSimulation_struct SSAL_RealisationSimulation;

/*macro functions*/


/*converts variable name list to indexes*/
#define SSAL_VARS2INDS(sim,model,varInd) {                  \
    int i,j;                                                \
    if ((sim)->Nvar == (model)->N)                          \
    {                                                       \
        for (j=0;j<(model)->N;j++)                          \
        {                                                   \
            varInd[j] = j;                                  \
        }                                                   \
        for (i=0;i<(sim)->Nvar;i++)                         \
        {                                                                \
            strncpy((sim)->var[i],(model)->names[i],SSAL_MAX_NAME_SIZE); \
        }                                                                \
    }                                                       \
    for (i=0;i<(sim)->Nvar;i++)                             \
    {                                                       \
        for (j=0;j<(model)->N;j++)                          \
        {                                                   \
            if (!strcmp((sim)->var[i],(model)->names[j]))   \
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
void SSAL_InitRNGS(unsigned int *, int);

void SSAL_HandleError(int , char *,int, unsigned char, unsigned char, char *); 

int SSAL_Initialise(int,char **); 

/* generic wrapper functions*/
SSAL_Model SSAL_ImportLSBML(const char * );
SSAL_Simulation SSAL_CreateRealisationsSim(SSAL_Model *,int, char **, int, int, 
                                SSAL_real_t *, SSAL_real_t * );

/* simulation and operations*/
int SSAL_Simulate(SSAL_Simulation *, SSAL_AlgorithmType, const char *);
/*output*/
int SSAL_WriteSimulation(FILE *,SSAL_Simulation);
int SSAL_WriteRealisationsSim(FILE *, SSAL_RealisationSimulation *);


/*==================================================================*/
/*         Chemical Reaction Network API                            */
/*==================================================================*/

/*object creation*/
SSAL_Model SSAL_CreateChemicalReactionNetwork(char **, int ,int, 
                            SSAL_real_t * restrict , SSAL_real_t * restrict , SSAL_real_t * restrict );

/* Chemical Reaction Network structures and functions*/
int SSAL_SimulateCRN(SSAL_Simulation *, SSAL_AlgorithmType, const char *);
int SSAL_SimulateCRNRealisations(SSAL_RealisationSimulation *, 
            SSAL_ChemicalReactionNetwork *, SSAL_AlgorithmType,  int, char **);
int SSAL_WriteChemicalReactionNetwork(FILE *,SSAL_ChemicalReactionNetwork);

/*==================================================================*/
/*         Stochastic Differential Equation API                     */
/*==================================================================*/
/* SDE structures and functions*/
SSAL_Model SSAL_CreateStochasticDifferentialEquation(char ** , int, int ,
          void (*)(SSAL_real_t *, uint32_t,SSAL_real_t*, uint32_t, SSAL_real_t,SSAL_real_t*),
          void (*)(SSAL_real_t *, uint32_t,SSAL_real_t*,uint32_t, SSAL_real_t,SSAL_real_t*), SSAL_real_t *);

int SSAL_SimulateSDE(SSAL_Simulation *, SSAL_AlgorithmType, const char *);
int SSAL_SimulateSDERealisations(SSAL_RealisationSimulation *, 
            SSAL_StochasticDifferentialEquation *, SSAL_AlgorithmType,  int, char **);

/*utilities*/
char** SSAL_UtilTokeniseArgs(int *,const char *);
#endif
