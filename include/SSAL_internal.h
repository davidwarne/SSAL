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
 * @file SSAL_internal.h
 * @brief Computational backend for SSAL
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Science and Engineering Faculty
 * @author Queensland University of Technology
 *
 */

#ifndef __SSAL_INTERNAL_H_
#define __SSAL_INTERNAL_H_

#ifdef __SERIAL__
#include "ESSA_sequential.h"
#include "util_sequential.h"
#include "ASSA_sequential.h"
#include "ODE_sequential.h"
#endif /* __SERIAL__ */

#ifdef __PARALLEL__
#include "ESSA_parallel.h"
#include "ASSA_parallal.h"
#include "util_parallel.h"
#endif /*__PARALLEL__*/

#endif /* __SSAL_INTERNAL_H_ */
