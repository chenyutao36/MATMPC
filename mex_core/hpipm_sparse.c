/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is mxFree software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, giaf (at) dtu.dk                                                       *
*                                                                                                 *
**************************************************************************************************/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"
#include "hpipm_d_ocp_qp_ipm.h"
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
		
	// get data 
    
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *C = mxGetPr( mxGetField(prhs[0], 0, "Cgx") );
    double *D = mxGetPr( mxGetField(prhs[0], 0, "Cgu") );
    double *CN = mxGetPr( mxGetField(prhs[0], 0, "CgN") );
    double *q = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *r = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *b = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *lg = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *ug = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *lbu = mxGetPr( mxGetField(prhs[0], 0, "lb_du") );
    double *ubu = mxGetPr( mxGetField(prhs[0], 0, "ub_du") );
    double *lbx = mxGetPr( mxGetField(prhs[0], 0, "lb_dx") );
    double *ubx = mxGetPr( mxGetField(prhs[0], 0, "ub_dx") );
    
    double *x = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *pi = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u_new") );
    double *mu_x = mxGetPr( mxGetField(prhs[0], 0, "mu_x_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int ng = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ngN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") ); 
    int nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
    
    double mu0 = mxGetScalar( mxGetField(prhs[0], 0, "mu0") ); 
    int max_qp_it = mxGetScalar( mxGetField(prhs[0], 0, "max_qp_it") );
    int pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "pred_corr") );
    int cond_pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "cond_pred_corr") );
    int solver_mode = mxGetScalar( mxGetField(prhs[0], 0, "solver_mode") );
               	
	int ii, jj;
    int idx;
    
    // number of states for each stage
	int nx_v[N+1];
	for(ii=0; ii<=N; ii++)
		nx_v[ii] = nx;

    // number of controls for each stage
	int nu_v[N+1];
	for(ii=0; ii<N; ii++)
		nu_v[ii] = nu;
	nu_v[N] = 0;

    // number of bounds for each stage
	int nb_v[N+1];
    nb_v[0] = nu+nx;
	for(ii=1; ii<N; ii++)
		nb_v[ii] = nu+nbx;
    nb_v[N] = nbx;
    
    // number of control bounds for each stage
	int nbu_v[N+1];
	for(ii=0; ii<N; ii++)
        nbu_v[ii] = nu;
	nbu_v[N] = 0;

    // number of state bounds for each stage
	int nbx_v[N+1];
    nbx_v[0] = nx;
	for(ii=1; ii<=N; ii++)
        nbx_v[ii] = nbx;
		
    // number of general bounds for each stage
	int ng_v[N+1];
	for(ii=0; ii<N; ii++)
		ng_v[ii] = ng;
	ng_v[N] = ngN;

	int ns_v[N+1];
	for(ii=0; ii<=N; ii++)
		ns_v[ii] = 0;
    
    int nsbx_v[N+1];
    int nsbu_v[N+1];
    int nsbg_v[N+1];
	for(ii=0; ii<=N; ii++){
		nsbx_v[ii] = 0;
        nsbu_v[ii] = 0;
        nsbg_v[ii] = 0;
    }
    

	int *hidxb[N+1];
    for(ii=0; ii<=N; ii++)
		hidxb[ii] = (int *)mxCalloc(nb_v[ii],sizeof(int));
    
    for(jj=0; jj<nbu_v[0]; jj++)
		hidxb[0][jj] = jj;
    for(jj=0; jj<nbx_v[0]; jj++)
		hidxb[0][nbu_v[0]+jj] = nbu_v[0]+jj;
	for(ii=1; ii<=N; ii++){		
		for(jj=0; jj<nbu_v[ii]; jj++)
			hidxb[ii][jj] = jj;
        
        for(jj=0; jj<nbx_v[ii]; jj++){
            idx = (int)nbx_idx[jj]-1;
			hidxb[ii][nbu_v[ii]+jj] = nbu_v[ii]+idx;
        }
	}	

	double *hA[N];
	double *hB[N];
	double *hb[N];
	double *hQ[N+1];
	double *hS[N];
	double *hR[N];
	double *hq[N+1];
	double *hr[N];
	double *hlb[N+1];
	double *hub[N+1];
	double *hC[N+1];
	double *hD[N];
	double *hlg[N+1];
	double *hug[N+1];
	double *hx[N+1];
	double *hu[N+1];
	double *hpi[N];
	double *hlam_lb[N+1];
	double *hlam_ub[N+1];
	double *hlam_lg[N+1];
	double *hlam_ug[N+1];

    for(ii=0; ii<N; ii++)
		hA[ii] = A+ii*nx*nx;

	for(ii=0; ii<N; ii++)
		hB[ii] = B+ii*nx*nu;
    
	for(ii=0; ii<N; ii++)
		hb[ii] = b+ii*nx;

	for(ii=0; ii<=N; ii++)
		hQ[ii] = Q+ii*nx*nx;
        
	for(ii=0; ii<N; ii++)
		hS[ii] = S+ii*nx*nu;

    for(ii=0; ii<N; ii++)
		hR[ii] = R+ii*nu*nu;

	for(ii=0; ii<=N; ii++)
		hq[ii] = q+ii*nx;
	
	for(ii=0; ii<N; ii++)
		hr[ii] = r+ii*nu;

    hlb[0] = (double *)mxCalloc(nb_v[0],sizeof(double));
    hub[0] = (double *)mxCalloc(nb_v[0],sizeof(double));
    for(ii=0; ii<nu; ii++){
		hlb[0][ii] = lbu[ii];
        hub[0][ii] = ubu[ii];
    }
    for(ii=0; ii<nx; ii++){
		hlb[0][nu+ii] = ds0[ii];
        hub[0][nu+ii] = ds0[ii];
    }
	for(ii=1; ii<=N; ii++){
        hlb[ii] = (double *)mxCalloc(nb_v[ii],sizeof(double));
        hub[ii] = (double *)mxCalloc(nb_v[ii],sizeof(double));
        for(jj=0;jj<nbu_v[ii];jj++){
            hlb[ii][jj] = lbu[ii*nu+jj];
            hub[ii][jj] = ubu[ii*nu+jj];
        }
        for(jj=0;jj<nbx_v[ii];jj++){
            hlb[ii][nbu_v[ii]+jj] = lbx[(ii-1)*nbx+jj];
            hub[ii][nbu_v[ii]+jj] = ubx[(ii-1)*nbx+jj];
        }
    }

	for(ii=0; ii<N; ii++)
		hC[ii] = C+ii*ng*nx;
	hC[N] = CN;

	for(ii=0; ii<N; ii++)
		hD[ii] = D+ii*ng*nu;

	for(ii=0; ii<=N; ii++)
		hlg[ii] = lg+ii*ng;

	for(ii=0; ii<=N; ii++)
		hug[ii] = ug+ii*ng;
	
	for(ii=0; ii<=N; ii++)
		hx[ii] = x+ii*nx;

	for(ii=0; ii<N; ii++)
		hu[ii] = u+ii*nu;
	
	for(ii=0; ii<N; ii++)
		hpi[ii] = pi+(ii+1)*nx;
    
    for(ii=0; ii<=N; ii++){
		hlam_lb[ii] = (double *)mxCalloc(nb_v[ii],sizeof(double));
		hlam_ub[ii] = (double *)mxCalloc(nb_v[ii],sizeof(double));
		hlam_lg[ii] = (double *)mxCalloc(ng_v[ii],sizeof(double));
        hlam_ug[ii] = (double *)mxCalloc(ng_v[ii],sizeof(double));
    }
    
	// qp dim
	int dim_size = d_memsize_ocp_qp_dim(N);
    void *dim_mem = mxCalloc(dim_size,1);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);
	d_cvt_int_to_ocp_qp_dim(N, nx_v, nu_v, nbx_v, nbu_v, ng_v, nsbx_v, nsbu_v, nsbg_v, &dim);
    
	// qp
	int qp_size = d_memsize_ocp_qp(&dim);
    void *qp_mem = mxCalloc(qp_size,1);
	struct d_ocp_qp qp;
	d_create_ocp_qp(&dim, &qp, qp_mem);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &qp);


	// qp sol
	int qp_sol_size = d_memsize_ocp_qp_sol(&dim);
    void *qp_sol_mem = mxCalloc(qp_sol_size,1);
	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(&dim, &qp_sol, qp_sol_mem);


	// ipm arg
	int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(&dim);
    void *ipm_arg_mem = mxCalloc(ipm_arg_size,1);
	struct d_ocp_qp_ipm_arg arg;
	d_create_ocp_qp_ipm_arg(&dim, &arg, ipm_arg_mem);
    
    // select the mode
    enum hpipm_mode mode;
    switch (solver_mode)
    {
        case 0: 
            mode = SPEED_ABS; 
            break;
        case 1:
            mode = SPEED; 
            break;
        case 2:
            mode = BALANCE; 
            break;
        case 3:
            mode = ROBUST; 
            break;
        default:
            mode = SPEED; 
    }   
	d_set_default_ocp_qp_ipm_arg(mode, &arg);

	arg.alpha_min = 1e-8;
	arg.res_g_max = 1e-4;
	arg.res_b_max = 1e-6;
	arg.res_d_max = 1e-6;
	arg.res_m_max = 1e-6;
	arg.mu0 = mu0;
	arg.iter_max = max_qp_it;
	arg.stat_max = max_qp_it;
	arg.pred_corr = pred_corr;
	arg.cond_pred_corr = cond_pred_corr;


	// ipm
	int ipm_size = d_memsize_ocp_qp_ipm(&dim, &arg);
    void *ipm_mem = mxCalloc(ipm_size,1);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim, &arg, &workspace, ipm_mem);

	// call solver
	int hpipm_return = d_solve_ocp_qp_ipm(&qp, &qp_sol, &arg, &workspace);
        
    // print stats
    int err = 0;
    if (workspace.qp_res[0]>arg.res_g_max){
        mexPrintf("res_g:%5.3e   res_g_max:%5.3e\n", workspace.qp_res[0], arg.res_g_max);
        err++;
    }
    if (workspace.qp_res[1]>arg.res_b_max){
        mexPrintf("res_b:%5.3e   res_b_max:%5.3e\n", workspace.qp_res[1], arg.res_b_max);
        err++;
    }
    if (workspace.qp_res[2]>arg.res_d_max){
        mexPrintf("res_g:%5.3e   res_g_max:%5.3e\n", workspace.qp_res[2], arg.res_d_max);
        err++;
    }
    if (workspace.qp_res[3]>arg.res_m_max){
        mexPrintf("res_g:%5.3e   res_g_max:%5.3e\n", workspace.qp_res[3], arg.res_m_max);
        err++;
    }   
    if (err > 0)
        mexErrMsgTxt("QP solver does not converge!");
    
    if (hpipm_return==1)
        mexErrMsgTxt("QP solver reaches maximum number of iterations!");

	// convert back solution
	d_cvt_ocp_qp_sol_to_colmaj(&qp_sol, hu, hx, NULL, NULL, hpi, hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);
    
    // extrac multipliers
    for(jj=0;jj<nx;jj++)
        pi[jj] = hlam_ub[0][nu+jj] - hlam_lb[0][nu+jj];
    
    for(ii=0;ii<N;ii++){       
        for(jj=0;jj<nu;jj++)
            mu_u[ii*nu+jj] = hlam_ub[ii][jj] - hlam_lb[ii][jj];
                
        for(jj=0;jj<ng;jj++)
            mu[ii*ng+jj] = hlam_ug[ii][jj] - hlam_lg[ii][jj];
    }
    for(jj=0;jj<ngN;jj++)
        mu[N*ng+jj] = hlam_ug[N][jj] - hlam_lg[N][jj];
    
    for(ii=0;ii<N;ii++){  
        for(jj=0;jj<nbx;jj++)
            mu_x[ii*nbx+jj] = hlam_ub[ii+1][nbu_v[ii+1]+jj] - hlam_lb[ii+1][nbu_v[ii+1]+jj];
    }       
                     
    // Free memory
	for(ii=0;ii<=N;ii++){
        mxFree(hidxb[ii]);
        mxFree(hlb[ii]);
        mxFree(hub[ii]);
        mxFree(hlam_lb[ii]);
        mxFree(hlam_ub[ii]);
        mxFree(hlam_lg[ii]);
        mxFree(hlam_ug[ii]);
    }
    
	mxFree(dim_mem);
	mxFree(qp_mem);
	mxFree(qp_sol_mem);
	mxFree(ipm_arg_mem);
	mxFree(ipm_mem);

}

