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
 */

#include "SSAL.h"

/**
 * @brief Handle error reporting
 * @param errCode The Error code 
 * @param funcName The Function namee in which the error occurred
 * @param lineNum The line number at which the error occurred
 * @param fatal Flag indicating the error was fatal and the progam should be 
 *              aborted
 * @param warning Indicates the error is to be interpreted as a warning
 * @param msg A custom message which can be appended to the default
 */
void 
SSAL_HandleError(int errCode, char * funcName, int lineNum, unsigned char fatal, 
                 unsigned char warning, char *msg)
{
    char * errName;
    char * defaultMsg;
    char * errType;
    switch (errCode)
    {
        case SSAL_UNKNOWN_OPTION_ERROR:
            errName = "SSAL_UNKNOWN_OPTION_ERROR";
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
int 
SSAL_Initialise(int argc,char **argv)
{
    int i;
    unsigned char infoflag;
    unsigned int seed;
    seed = 1337;
    infoflag = 0;
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
    /*parse args*/
    for (;i<argc;i++)
    {
        if (!strcmp("--seed",argv[i]))
        {
            seed = (unsigned int)atoi(argv[++i]);          
        }
        else if (!strcmp("--info",argv[i]))
        {
            infoflag = 1;
        }
        else
        {
            SSAL_HandleError(SSAL_UNKNOWN_OPTION_ERROR,"SSAL_Initialise",
                             __LINE__,0,1,argv[i]);
        }
    }

    /*initialise RNGS */
    #ifdef __SERIAL__
       suarngs(seed,NULL,NULL);
    #endif
    if (infoflag)
    {
        fprintf(stderr,"RNG Seed: %d\n",seed);
        fprintf(stderr,"RNG RAND_MAX: %d\n",RAND_MAX);
    }
    return SSAL_SUCCESS;
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
 * an M by N matrix using row-major storage. nu[j][i] is the quantity of 
 * chemical species i required for reaction j.
 *
 * @note c must an M length vector, where c[j] is the kinetic rate of 
 * reaction j.
 *
 * @note Internally this function creates new copies of the stochiometric 
 * coefficient matrices. This ensures that the resulting allocated memory is 
 * valid, however it does not prevent invalid data in the matrices if the input 
 * pointers are not valid.
 */
SSAL_CRN 
SSAL_CreateChemicalReactionNetwork(char ** names, int M,int N, 
                                   SSAL_real_t * restrict nu_minus, 
                                   SSAL_real_t * restrict nu_plus, 
                                   SSAL_real_t * restrict c)
{
    char * funcname;
    char * nameArray;
    size_t NxM;
    SSAL_CRN newCRN;
   
   funcname = "SSAL_CreateChemicalReactionNetwork";
    
    /*allocate memory*/
    newCRN.N = (uint32_t)N;
    newCRN.M = (uint32_t)M;
    NxM = newCRN.N*newCRN.M;
    newCRN.nvar = N;
    if ((newCRN.vars = (int *)malloc(newCRN.N*sizeof(int))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }

    if ((nameArray = (char *)malloc(newCRN.N*SSAL_MAX_NAME_SIZE*sizeof(char)))
        == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    if((newCRN.names = (char **)malloc(newCRN.N*sizeof(char*))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    for (int i=0;i<newCRN.N;i++)
    {
        newCRN.names[i] = nameArray + i*SSAL_MAX_NAME_SIZE;
    }
    
    if((newCRN.nu_minus = (SSAL_real_t *)malloc(NxM*sizeof(SSAL_real_t))) 
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN.nu_plus = (SSAL_real_t *)malloc(NxM*sizeof(SSAL_real_t))) 
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    if((newCRN.nu = (SSAL_real_t *)malloc(NxM*sizeof(SSAL_real_t))) 
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }

    if((newCRN.c = (SSAL_real_t *)malloc(newCRN.M*sizeof(SSAL_real_t))) 
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    /*copy data to valid arrays*/
    for (int i=0;i<NxM;i++)
    {
        newCRN.nu_minus[i] = nu_minus[i];
    }
    for (int i=0;i<NxM;i++)
    {
        newCRN.nu_plus[i] = nu_plus[i];
    }
    for (int i=0;i<NxM;i++)
    {
        newCRN.nu[i] = nu_plus[i] - nu_minus[i];
    }

    for (int i=0;i<newCRN.M;i++)
    {
        newCRN.c[i] = c[i];
    }
    for (int i=0;i<newCRN.N;i++)
    {
        for (int j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newCRN.names[i][j] = names[i][j];
        }
        newCRN.names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    for (int i=0;i<newCRN.N;i++)
    {
        newCRN.vars[i] = i;
    }

    /*all done, return the new model*/
    return newCRN;
}

/**
 * @brief Safely creates a multivariate stochastic differential equation
 * @detail checks inputs are valid and creates the SDE
 *
 * @param M the number of parameters  
 * @param N  the number of equations  
 * @param mu function pointer to drift function
 * @param sigma function pointer diffusion function
 * @param p model parameter vector
 *
 * @returns A Stochastic Differential Equation Structure
 */
SSAL_SDE 
SSAL_CreateStochasticDifferentialEquation(char ** names, int M, int N, 
                                          void (*mu)(SSAL_real_t *, uint32_t, 
                                                    SSAL_real_t*, uint32_t,
                                                    SSAL_real_t,SSAL_real_t*), 
                                          void (*sigma)(SSAL_real_t *, uint32_t, 
                                                      SSAL_real_t*, uint32_t, 
                                                      SSAL_real_t,SSAL_real_t*), 
                                          SSAL_real_t *p)
{
    char * funcname;
    char * nameArray;
    SSAL_SDE newSDE;
   
   funcname = "SSAL_CreateStochasticDifferentialEquation";
   
    /*allocate memory*/
    newSDE.N = (uint32_t)N;
    newSDE.M = (uint32_t)M;
    
    if ((nameArray = (char *)malloc(newSDE.N*SSAL_MAX_NAME_SIZE*sizeof(char))) 
        == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newSDE.names = (char **)malloc(newSDE.N*sizeof(char*)))== NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    for (int i=0;i<newSDE.N;i++)
    {
        newSDE.names[i] = nameArray + i*SSAL_MAX_NAME_SIZE;
    }

    if ((newSDE.p = (SSAL_real_t *)malloc(newSDE.M*sizeof(SSAL_real_t))) 
        == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }

    for (int i=0;i<newSDE.M;i++)
    {
       newSDE.p[i] = p[i]; 
    }
   
    newSDE.mu = mu;
    newSDE.sigma = sigma;
    
    for (int i=0;i<newSDE.N;i++)
    {
        for (int j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newSDE.names[i][j] = names[i][j];
        }
        newSDE.names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    /*all done, return the new model*/
    return newSDE;
}

/**
 * @brief write chemical reation network too file
 * @param stream the output stream
 * @param data a chemical reaction network
 *
 */
int 
SSAL_WriteChemicalReactionNetwork(FILE * stream ,
                                  SSAL_CRN data)
{
    int i,j;

    if (stream == NULL)
    {
        return SSAL_IO_ERROR;
    }
    fprintf(stream,"{\"version\" : \"1\", \"level\" : \"3\",");
    fprintf(stream,"\"model\" : { \"id\" : \"crn\",\"name\" : \"crn\",\"substanceUnits\" : \"item\", \"timeUnits\" : \"seconds\",\"extentUnits\" : \"item\",");
    fprintf(stream,"\"compartments\" : [");
    fprintf(stream,"{\"id\" : \"space\",\"spatialDimensions\" : 0, \"size\" : 0}");
    fprintf(stream,"],");
    fprintf(stream,"\"species\" : [");
        fprintf(stream,"{\"id\" : \"%s\", \"compartment\" : \"space\",\"initialAmount\" : %d, \"hasOnlySubstanceUnits\" : true}",data.names[0],(int)data.X0[0]);
    for (i=1;i<data.N;i++)
    {
        fprintf(stream,",{\"id\" : \"%s\", \"compartment\" : \"space\",\"initialAmount\" : %d, \"hasOnlySubstanceUnits\" : true}",data.names[i],(int)data.X0[i]);
    }
    fprintf(stream,"],");
    fprintf(stream,"\"parameters\" : [");
        fprintf(stream,"{\"id\" : \"c%d\",\"value\" : %f}",0,data.c[0]);
    for (j=1;j<data.M;j++)
    {
        fprintf(stream,",{\"id\" : \"c%d\",\"value\" : %f}",j,data.c[j]);
    }
    fprintf(stream,"],");
    fprintf(stream,"\"reactions\" : [");
        fprintf(stream,"{\"id\" : \"R%d\",\"reversible\" : false, \"reactants\" : [",0);
        int first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_minus[i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
                fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_minus[i]);
            }
        }
        fprintf(stream,"],\"products\" : [");
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_plus[i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");

            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_plus[i]);
            }
        }
        fprintf(stream,"]}");
    for (j=1;j<data.M;j++)
    {
        fprintf(stream,",{\"id\" : \"R%d\",\"reversible\" : false, \"reactants\" : [",j);
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_minus[j*data.N + i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_minus[j*data.N + i]);
            }
        }
        fprintf(stream,"],\"products\" : [");
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_plus[j*data.N + i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_plus[j*data.N + i]);
            }
        }
        fprintf(stream,"]}");

    }
    fprintf(stream,"]");
    fprintf(stream,"}");
    fprintf(stream,"}\n");
    return SSAL_SUCCESS; 
}


/**
 * @brief imports chemical reaction network from an LSBML file.
 * @param filename the name of the lsbml file
 * @return A Chemical Reaction Network
 */
SSAL_CRN 
SSAL_ImportLSBML(const char * filename )
{
    cJSON *root, *model, *parameters, *reactions, *species;
    FILE *fp;
    long int size;
    char * lsbml_json;
    SSAL_CRN CRN;
    char * funcname;
    funcname = "SSAL_ImportLSBML";

    /*open the file*/
    if (!(fp = fopen(filename,"r")))
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    /*get file size in bytes*/
    fseek(fp,0,SEEK_END);
    size = ftell(fp);
    rewind(fp);
    /*allocate buffer*/
    if((lsbml_json = (char *)malloc(size*sizeof(char))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    /*read contents*/
    fread((void *)lsbml_json,sizeof(char),size,fp);

    /*parse the file*/
    root = cJSON_Parse(lsbml_json);
    model = cJSON_GetObjectItem(root,"model");
    parameters = cJSON_GetObjectItem(model,"parameters");
    reactions = cJSON_GetObjectItem(model,"reactions");
    species = cJSON_GetObjectItem(model,"species");

    CRN.M = cJSON_GetArraySize(reactions);
    CRN.N = cJSON_GetArraySize(species);  

    /*by default we set nvar and vars such that all species are tracked*/
    CRN.nvar = CRN.N;
    if ((CRN.vars = (int *)malloc(CRN.N*sizeof(int))) == NULL) 
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    for (int i=0;i<CRN.N;i++)
    {
        CRN.vars[i] = i;
    }
    
    /*allocate memory for CRN components*/
    if((CRN.names = (char **)malloc((CRN.N)*sizeof(char*))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    for (int i=0;i<CRN.N;i++)
    {
        if ((CRN.names[i] = (char *)malloc(SSAL_MAX_NAME_SIZE*sizeof(char))) 
             == NULL)
        {
            SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
        }
    }
    if((CRN.X0 = (SSAL_real_t*) malloc((CRN.N)*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((CRN.c = (SSAL_real_t*) malloc((CRN.M)*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((CRN.nu_minus = 
       (SSAL_real_t*)malloc((CRN.N)*(CRN.M)*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((CRN.nu_plus = (SSAL_real_t*)malloc((CRN.N)*(CRN.M)*sizeof(SSAL_real_t)))
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((CRN.nu = (SSAL_real_t*)malloc((CRN.N)*(CRN.M)*sizeof(SSAL_real_t))) 
       == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    memset(CRN.nu_minus,0,(CRN.N)*(CRN.M)*sizeof(SSAL_real_t));
    memset(CRN.nu_plus,0,(CRN.N)*(CRN.M)*sizeof(SSAL_real_t));
    /*read species names*/
    for (int i=0;i<CRN.N;i++)
    {
        strncpy(CRN.names[i],
           cJSON_GetObjectItem(cJSON_GetArrayItem(species,i),"id")->valuestring,
           SSAL_MAX_NAME_SIZE);
        CRN.X0[i] = (SSAL_real_t)(cJSON_GetObjectItem(
                                                  cJSON_GetArrayItem(species,i),
                                                  "initialAmount")->valuedouble);
    }
    /*read rate parameters*/
    /**@note we assume that the reaction rates are the only parameters, 
     * and that they are listed in the same order as the matching reaction
     */
    for (int j=0;j<CRN.M;j++)
    {
        CRN.c[j] = (SSAL_real_t)(cJSON_GetObjectItem(
                                               cJSON_GetArrayItem(parameters,j),
                                               "value")->valuedouble);
    }

    /* the reactions*/
    for (int j=0;j<CRN.M;j++)
    {
       cJSON *react, *prod;
       int nReact, nProd;
       char temp[SSAL_MAX_NAME_SIZE];
       react = cJSON_GetObjectItem(cJSON_GetArrayItem(reactions,j),"reactants");
       prod = cJSON_GetObjectItem(cJSON_GetArrayItem(reactions,j),"products");

        nReact = cJSON_GetArraySize(react);
        nProd = cJSON_GetArraySize(prod);
        for (int k=0;k<nReact;k++)
        {
            for (int i=0;i<CRN.N;i++)
            {
                if (!strncmp(CRN.names[i],
                             cJSON_GetObjectItem(cJSON_GetArrayItem(react,k),
                                                 "species")->valuestring,
                             SSAL_MAX_NAME_SIZE))
                {
                    CRN.nu_minus[j*(CRN.N)+i] = 
                             (SSAL_real_t)cJSON_GetObjectItem(
                                                  cJSON_GetArrayItem(react,k),
                                                  "stoichiometry")->valuedouble;
                }
            }
        }
        
        for (int k=0;k<nProd;k++)
        {
            for (int i=0;i<CRN.N;i++)
            {
                if (!strncmp(CRN.names[i],
                            cJSON_GetObjectItem(cJSON_GetArrayItem(prod,k),
                                                "species")->valuestring,
                            SSAL_MAX_NAME_SIZE))
                {
                    CRN.nu_plus[j*(CRN.N)+i] = 
                               (SSAL_real_t)cJSON_GetObjectItem(
                                                  cJSON_GetArrayItem(prod,k),
                                                  "stoichiometry")->valuedouble;
                }
            }
        }

    }

    for (int i=0;i<CRN.N*CRN.M;i++)
    {
        CRN.nu[i] = CRN.nu_plus[i] - CRN.nu_minus[i];
    }

    cJSON_Delete(root);
    return CRN;
}


/**
 * @brief utilitiy function which breaks a char array into a array of char arrays
 * @param argc arg count
 * @param args input arg string 
 * @return argv a array of pointers to a char arrays
 */
char** SSAL_UtilTokeniseArgs(int *argc,const char * args)
{
    char *buf;
    char **argv;

    int numChars;
    int numArgs;
    int argsLength;
    int i,j,k;
    char cur_char;
    char prev_char;
    char * funcname;
    funcname = "SSAL_UtilTokeniseArgs";
    cur_char = args[0];
    prev_char = ' ';
    numChars = 0;
    numArgs = 0;
    i = 0;
    if (args[0] == '\0')
    {
        *argc = 0;
        return NULL;
    }
    /*first pass to get sizes. */
    while (cur_char != '\0' && i != SSAL_MAX_BUFFER_SIZE)
    {
        numChars += (cur_char != ' ');
        numArgs += (cur_char != ' ' && prev_char == ' ');
        i++;
        prev_char = cur_char;
        cur_char = args[i];
    }


    argsLength = i; 
    if((buf = (char *)malloc((numChars+numArgs)*sizeof(char)))==NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if ((argv = (char **)malloc(numArgs*sizeof(char*))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    j = 0;
    k = 0;
    /*second pass to collect tokenised data*/
    for (i=0;i<argsLength;i++)
    {
        /*starting a token*/
        if (args[i] == ' ' && args[i+1] != ' ')
        {
            j++;
            argv[k] = buf + j;
            k++;
        } /*ending a token*/
        else if (i==0 && args[i] != ' ')
        {
            buf[j] = args[i];
            argv[k] = buf + j;
            j++;
            k++;
        }
        else if ( args[i] != ' ')
        {
            buf[j] = args[i];
            j++;
        }
    }

    /*append null characters*/
    for (k=1;k<numArgs;k++)
    {
        argv[k][-1] = '\0';
    }
    buf[numChars+numArgs - 1] = '\0';
    *argc = numArgs;
    return argv;
}

