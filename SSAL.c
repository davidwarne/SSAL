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
 * @file SSAL.c
 * @brief A fast implementation of standard Stochastic Simulation algorithms.
 * @details This library provides optimised sequential and parallel exact and 
 * approximate stochastic simulation algorithms. The implementations provide an
 * easy to use application program interface (API) using standard C structures.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Science and Engineering Faculty
 * @author Queensland University of Technology
 *
 * @date 10 Sep 2015
 */

#include "SSAL.h"
#include "SSAL_default_callbacks.h"

SSAL_Options SSAL_GLOBAL_Options;
SSAL_Algorithm SSAL_ALG_LUT[SSAL_NUM_ALGS];
unsigned char SSAL_ALG_ENABLED[SSAL_NUM_ALGS];
/**
 * @brief Handle error reporting
 * @param errCode The Error code 
 * @param funcName The Function namee in which the error occurred
 * @param lineNum The line number at which the error occurred
 * @param fatal Flag indicating the error was fatal and the progam should be aborted
 * @param warning Indicates the error is to be interpreted as a warning
 * @param msg A custom message which can be appended to the default
 */
void SSAL_HandleError(int errCode, char * funcName, int lineNum, unsigned char fatal, 
                    unsigned char warning, char *msg)
{
    char * errName;
    char * defaultMsg;
    char * errType;
    switch (errCode)
    {
        case SSAL_UNKNOWN_OPTION:
            errName = "SSAL_UNKNOWN_OPTION";
            defaultMsg = "A configuration option was not recognised.";
            break;
        default:
            errName = "SSAL_UNKNOWN_ERROR";
            defaultMsg = "No idea what happened then...";
            break;
    }

    if (warning)
    {
        errType = "ERROR";   
    }
    else
    {
        errType = "WARNING";
    }
        
    if (msg != NULL)
    {
        fprintf(stderr,"[%s %s]: %s\n %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,msg,errCode,lineNum);    
    }
    else
    {
        fprintf(stderr,"[%s %s]: %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,errCode,lineNum);    
    }
    

    if (fatal)
    {
        fprintf(stderr,"Fatal Error Aborting!\n");
        exit(1);
    }
    return;
}



/**
 * @brief Initilialises the SSAL library
 * @detail Initialised the computational backend based on compilation
 * options applied. Some user customisations are also available.
 *
 * @param argc the number of input args
 * @param argv an array of args
 * @todo Will add more options here as needed.
 * @retVal SSAL_SUCCESS if initialisation was successful
 */
int SSAL_Inititalise(int argc,char **argv)
{
    int i;
    /*search for config options*/
    for (i=1;i<argc;i++)
    {
        if (!strcmp("--SSALOPT",argv[i]))
        {
            i++;

            break;
        }
    }

    /*set defaults*/
    SSAL_GLOBAL_Options.algOptions = SSAL_ALG_DEFAULT;
    
    /*parse args*/
    for (;i<argc;i++)
    {
        if (!strcmp("--sequential-Only",argv[i]))
        {
            SSAL_GLOBAL_Options.algOptions = SSAL_ALG_SEQUENTIAL;
        }
        else
        {
            SSAL_HandleError(SSAL_UNKNOWN_OPTION,"SSAL_Initialise",__LINE__,0,1,argv[i]);
        }
    }
    /*register algorithms*/
    /*set all callback to NULL*/
    for (i=0;i<SSAL_NUM_ALGS;i++)
    {
        SSAL_REGISTER_ALG(i,NULL,NULL,NULL,NULL)
    }
    /*register callbacks*/
    SSAL_REGISTER_ALG(SSAL_ESSA_GILLESPIE_SEQUENTIAL,&checkANYCRNSEQ,&exec_essa_gil_seq,NULL,NULL)
    #ifdef __PARALLEL__
    #endif
}

/**
 * @brief Safely creates a chemical reaction network
 * @detail checks inputs are valid and creates the network
 *
 * @param N the number of chemical species
 * @param M the number of chemical reactions
 * @param nu_minus the reactant stochiometric coefficient matrix
 * @param nu_plus the product stochiometric coefficient matrix
 * @param c kinetic rate contstants
 *
 * @returns A Chemical Reaction Network Structure
 *
 * @note nu_minus and nu_minus must be N*M length vectors representing 
 * an M by N matrix using row-major storage. nu[j][i] is the quantity of chemical
 * species i required for reaction j.
 *
 * @note c must an M length vector, where c[j] is the kinetic rate of reaction j.
 *
 * @note Internally this function creates new copies of the stochiometric 
 * coefficient matrices. This ensures that the resulting allocated memory is valid, 
 * however it does not prevent invalid data in the matrices if the input pointers 
 * are not valid.
 */
SSAL_Model SSAL_CreateChemicalReactionNetwork(char ** names, int M,int N, 
                            float * restrict nu_minus, float * restrict nu_plus, float * restrict c)
{
    int i;
    char * nameArray;
    SSAL_Model newModel;
    SSAL_ChemicalReactionNetwork *newCRN;
   
    newModel.type = SSAL_CHEMICAL_REACTION_NETWORK;
    newCRN = (SSAL_ChemicalReactionNetwork*)malloc(sizeof(SSAL_ChemicalReactionNetwork));
    newModel.model = (void *)newCRN;
    /*allocate memory*/
    newCRN->N = (uint32_t)N;
    newCRN->M = (uint32_t)M;
    nameArray = (char *)malloc(newCRN->N*SSAL_MAX_NAME_SIZE*sizeof(char));
    newCRN->names = (char **)malloc(newCRN->N*sizeof(char*));
    for (i=0;i<newCRN->N;i++)
    {
        newCRN->names[i] = nameArray + i*SSAL_MAX_NAME_SIZE;
    }
    newCRN->nu_minus = (float *)malloc(newCRN->N*newCRN->M*sizeof(float));
    newCRN->nu_plus = (float *)malloc(newCRN->N*newCRN->M*sizeof(float));
    newCRN->c = (float *)malloc(newCRN->M*sizeof(float));
    
    /*copy data to valid arrays*/
    for (i=0;i<newCRN->N*newCRN->M;i++)
    {
        newCRN->nu_minus[i] = nu_minus[i];
    }
    for (i=0;i<newCRN->N*newCRN->M;i++)
    {
        newCRN->nu_plus[i] = nu_plus[i];
    }
    for (i=0;i<newCRN->M;i++)
    {
        newCRN->c[i] = c[i];
    }
    for (i=0;i<newCRN->N;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newCRN->names[i][j] = names[i][j];
            if (names[i][j] == '\0')
            {
                continue;
            }
        }
        newCRN->names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    /*all done, return the new model*/
    return newModel;
}

/**
 * @brief write chemical reation network too file
 * @param stream the output stream
 * @param data a chemical reaction network
 */
int SSAL_WriteChemicalReactionNetwork(FILE * stream ,SSAL_ChemicalReactionNetwork data)
{
    int i,j;
    if (stream == NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,"SSAL_WriteChemicalReactionNetwork",__LINE__,1,0,NULL);
    }

    /*print species data*/
    fprintf(stream,"Chemical Species:\n");
    for (i=0;i<data.N;i++)
    {
        fprintf(stream,"\tX[%d]: %s\n",i,data.names[i]);    
    }

    /*print equations*/
    fprintf(stream,"Chemical Reactions:\n");
    fprintf(stream,"\t         c \n");
    fprintf(stream,"\t nu^-*X --> nu^+*X\n");
    fprintf(stream,"\tnu^- =\n");
    for (j=0;j<data.M;j++)
    {
        fprintf(stream,"\t\t |");
        for (i=0;i<data.N;i++)
        {
            fprintf(stream," %f",data.nu_minus[j*data.N + i]);
        }
        fprintf(stream," |\n");
    }
    fprintf(stream,"\tnu^+ =\n");
    for (j=0;j<data.M;j++)
    {
        fprintf(stream,"\t\t |");
        for (i=0;i<data.N;i++)
        {
            fprintf(stream," %f",data.nu_minus[j*data.N + i]);
        }
        fprintf(stream," |\n");
    }
    fprintf(stream,"Kinetic Reaction Rates:\n");
    for (i=0;i<data.M;i++)
    {
        fprintf(stream,"\tc[%d] = %f\n",i,data.c[i]);
    }
    return SSAL_SUCCESS; 
}

/**
 * @brief Creates a realisation simulation structure
 * @detail pre-allocates memory in an efficient structure for the model and algorithm
 * @param model a pointer referencing the model realisations will be made of
 * @param N number of variables to observe
 * @param obs a pointer of model entities to observe
 * @param NR number of realisations
 * @param NT number of time points to record
 * @param T time points to record
 * @param initCond initial conditions 
 */
SSAL_Simulation SSAL_CreateRealisationsSim(SSAL_Model *model, int N,char **obs, int NR, int NT, 
                                float* T, float * initCond)
{
    int i;
    char * obsArray;
    SSAL_Simulation newSim;
    SSAL_RealisationSimulation *newRS;
    newSim.type = SSAL_REALISATIONS;
    newSim.model = model;
    newRS = (SSAL_RealisationSimulation*)malloc(sizeof(SSAL_RealisationSimulation));
    newSim.sim = (void*)newRS;
    /*allocate memory*/
    newRS->NR = NR;
    newRS->NT = NT;
    newRS->T = (float *)malloc(NT*sizeof(float));
    /**@todo extend to handle multiple scenarios */
    newRS->NIC = 1; 
    newRS->Nvar = N;
    obsArray = (char *)malloc(N*SSAL_MAX_NAME_SIZE*sizeof(char));
    newRS->var = (char **)malloc(N*sizeof(char*));
    for (i=0;i<newRS->Nvar;i++)
    {
        newRS->var[i] = obsArray + i*SSAL_MAX_NAME_SIZE;
    }

    /*copy data */
    for (i=0;i<newRS->NT;i++)
    {
        newRS->T[i] = T[i];
    }
    for (i=0;i<newRS->Nvar;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newRS->var[i][j] = obs[i][j];
            if (obs[i][j] == '\0');
            {
                continue;
            }
        }
        newRS->var[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }


    /*setting initial condition sizes are dependent on the model*/
    switch (model->type)
    {
        case SSAL_CHEMICAL_REACTION_NETWORK:
        {
            SSAL_ChemicalReactionNetwork *CRN;
            CRN = (SSAL_ChemicalReactionNetwork *)model->model;
            
            /*initial conditions will be the initial chemical species copy numbers*/
            newRS->initCond = (float *)malloc((newRS->Nvar)*(CRN->N)*sizeof(float));
            for (i=0; i<(newRS->Nvar)*(CRN->N); i++)
            {
                newRS->initCond[i] = initCond[i];
            }
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNKNOWN_TYPE,"SSAL_CreateRealisationsSimulation",
                            __LINE__,1,0,NULL);        
            break;
    }
    return newModel;
}


/**
 * @brief Runs a stochastic simulation
 *
 * @param sim pointer to simulation structure
 * @param alg the selected algorithm 
 * @param args algorithm specific parameters
 * 
 * @return return code
 */
int SSAL_Simulate(SSAL_Simulation *sim, SSAL_Algorithm alg, const char * args)
{

    if (alg < SSAL_NUM_ALGS && alg >= 0)
    {
        SSAL_INPUTCHECK(alg,sim);            
        SSAL_EXEC(alg,sim,nparams,params);            
    }
}


 
