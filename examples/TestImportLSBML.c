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
 * @file TestImportLSBML.c
 * @brief An example program using the SSAL Import function to load a
 * lightweight SBML file (i.e., a JSON version of s subset of SBML).
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Queensland University of Technology
 *
 */
#include<stdlib.h>
#include<stdio.h>
#include "SSAL.h"


int 
main(int argc , char ** argv)
{

    SSAL_CRN CRN;
    if (argc != 2)
    {
        fprintf(stderr,"Usage: %s <filename>\n",argv[0]);
    }

    /*init SSAL */
    SSAL_Initialise(argc,argv);

    CRN = SSAL_ImportLSBML(argv[1]);
    SSAL_WriteChemicalReactionNetwork(stderr,CRN);

    return 0;
}
