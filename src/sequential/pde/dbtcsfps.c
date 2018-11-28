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
#include "PDE_sequential.h"

/**
 * @brief Solve parabolic partial differential equation of the form,
 *        u_t = (D(u) u_x)_x + f(u)
 * on the spatial interval (0,L] over time (0,T]. Subject to zero flux initial
 * conditions.
 * @detail Uses the time-implicit backward in time/forward in space 
 * discretisation. This results in an equation of the form,
 *               u_i^{j+1} = G(u_i^{j+1}) + u_i^j
 * where 
 *   G(u_i^{j+1}) = delta_t/delta_x^2 [ 
 *                      D(u_{i+1/2}^{j+1})(u_{i+1}^{j+1} - u_i^{j+1}) 
 *                    - D(u_{i-1/2}^{j+1})(u_i^{j+1} - u_{i-1}^{j+1})]i
 *                  + delta_t*f(u_i^{j+1})
 *  and G(u_0^{j+1}) = G(u_{Nx-1}^{j+1}) = 0.
 *  The half steps are approximated with the arithmetic average 
 *  D(u_{i+1/2}^{j+1}) = (D(U_{i+1}^{j+1}) + D(u_i^{j+1}))/2
 *  Fixed point iteration is applied until the residuals are within a user 
 *  specified tolerance. If more than maxiters is performed the method returns 
 *  -1 to indicate the fixed-point iterations have not converged.
 *
 * @param L spatial domain [0,L]
 * @param T temporal domain [0,T]
 * @param nt number of observation times
 * @param Ti array of observation times
 * @param p model parameters
 * @param u0 array of size Nx storing u(x,0) at mesh points
 * @param D u dependent diffusion function
 * @param f u dependent source/sink function
 * @param Nx number of spatial nodes
 * @param Nt number of temporal nodes
 * @param tol tolerance on fixed-point iteration residual
 * @param maxiters convergence failed if maxiters fixed-point iterations fail to
 * obtain residuals within tol.
 * @param u_r PDE solution at observation times (nt*Nx) 
 * 
 * @returns 0 if successful solve, -1 if fixed-point iterations failed to 
 * converged
 * @todo extend to support other bounndary conditions
 */
int
dbtcsfps(double L, double T, int nt, double* restrict Ti, double * restrict p, 
         double * restrict u0, void (*D)(double *, int , double*, double *), 
         void (*f)(double *, int, double *, double *),int Nx, int Nt, double tol,
         int maxiters, double * restrict u_r)
{
    double dx; /*spatial step*/
    double dt; /*time-step*/
    double t; /*current time*/
    double *up; /*previous timestep*/
    double *u_buf; /*double buffer for fast swap in fixed-point iteration*/
    double *un; /*for fixed-point u_n = F(u_{n-1})*/
    double *unm1; /*u_{n-1}*/
    double *Du; /*diffusion */
    double *fu; /*source*/
    double D_bar; /*max diffusion for PDE stability test*/
    double u_test;

    /*generate mesh spacings*/
    dx = L/((double)(Nx - 1));
    dt = T/((double)(Nt - 1));
    /*@todo for now we assume max(D(u)) occurs at D(1)*/
    u_test = 1.0;
    (*D)(&u_test,1,p,&D_bar);
    /*check stability condition*/
//    if (dt*D_bar > dx*dx )
//    {
//        fprintf(stderr,"Warning: dt*D > dx*dx, so solution may not be stable. Either decrease dt or increase dx\n.");
//    }
    /*allocate memory for current and previous times steps*/
    up = (double*)malloc(Nx*sizeof(double));    
    u_buf = (double*)malloc(2*Nx*sizeof(double));    

    unm1 = u_buf;    
    un = u_buf + Nx;    
    fu = (double*)malloc(Nx*sizeof(double));    
    Du = (double*)malloc(Nx*sizeof(double));

    /*intialised*/
    memcpy(un,u0,Nx*sizeof(double));
    memcpy(unm1,u0,Nx*sizeof(double));
    memcpy(up,u0,Nx*sizeof(double));
    t = 0;

    /*solve for u at each observation time*/
    for (int ti=0;ti<nt;ti++)
    {
        /*time-step between observations times*/
        for(; t <= (Ti[ti] -dt); t+=dt)
        {
            /* fixed-point iterations till < tol*/
            double res;
            int iter;
            /*store copy of previous timestep*/                     
            memcpy(up,un,Nx*sizeof(double));
            
            res = 1.0;
            iter = 0;
            /*apply fixed point iterations*/
            while (res > tol && iter  < maxiters)
            {
                double *tmp_ptr;
                tmp_ptr = NULL;
                /*swap buffers*/
                tmp_ptr = un;
                un = unm1;
                unm1 = tmp_ptr;

                /*compute D and S at each node*/
                (*D)(unm1,Nx,p,Du);  
                (*f)(unm1,Nx,p,fu);
                un[0] = 0.0; 
                /* compute BTCS discretisation*/
                for (int i=1;i<Nx-1;i++)
                {
                    un[i] = 0.5*(dt/(dx*dx))*(
                                        (Du[i] + Du[i+1])*(unm1[i+1] - unm1[i]) 
                                      - (Du[i] + Du[i-1])*(unm1[i] - unm1[i-1])
                                      );
                }
                un[Nx-1] = 0.0;
                /* add source terms an previous solution*/
                for (int i=0;i<Nx;i++)
                {
                    un[i] += dt*fu[i] + up[i];
                }
                /*compute residual in l_2 norm for convergence check*/
                res = 0;
                for (int i=0;i<Nx;i++)
                {
                    res += (un[i] - unm1[i])*(un[i] - unm1[i]);
                }
                res = sqrt(res);
                iter++;
            }

            /*failed to converge*/
            if (res > tol)
            {
                fprintf(stderr,"Error: Failed to converge.\n");
                free(Du);
                free(fu);
                free(up);
                free(u_buf);
                return -1;
            }
        }
        /*store times step*/
        memcpy(u_r + ti*Nx,un,Nx*sizeof(double));
    }
    free(Du);
    free(fu);
    free(up);
    free(u_buf);
    return 0;
}
