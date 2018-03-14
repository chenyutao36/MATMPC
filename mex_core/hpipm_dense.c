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

#include "hpipm_d_dense_qp_dim.h"
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
#include "hpipm_d_dense_qp_ipm.h"
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
		
	// get data 
    
    double *Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc") );
    double *gc = mxGetPr( mxGetField(prhs[0], 0, "gc") );
    double *Cc = mxGetPr( mxGetField(prhs[0], 0, "Cc") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb_du") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub_du") );
    double *lg = mxGetPr( mxGetField(prhs[0], 0, "lcc") );
    double *ug = mxGetPr( mxGetField(prhs[0], 0, "ucc") );
    
    double *u = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
    
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int nc = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ncN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") );  
    
    double mu0 = mxGetScalar( mxGetField(prhs[0], 0, "mu0") ); 
    int max_qp_it = mxGetScalar( mxGetField(prhs[0], 0, "max_qp_it") );
    int pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "pred_corr") );
    int cond_pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "cond_pred_corr") );
        
    int nv = N*nu;
    int ne = 0;
    int nb = N*nu;
    int ng = N*nc+ncN;
    int ns = 0;
    
    plhs[0] = mxCreateDoubleMatrix(ng, 1, mxREAL);
    double *mu_vec = mxGetPr(plhs[0]);
    
    int i,j;
    int *idxb = (int *) mxMalloc(nb*sizeof(int));
    for (i=0;i<nb;i++)
        idxb[i] = i;
    
    // dense qp dim
    int dense_qp_dim_size = d_memsize_dense_qp_dim();
	void *dense_qp_dim_mem = mxCalloc(dense_qp_dim_size,1);

	struct d_dense_qp_dim qp_dim;
	d_create_dense_qp_dim(&qp_dim, dense_qp_dim_mem);

	d_cvt_int_to_dense_qp_dim(nv, ne, nb, ng, ns, &qp_dim);
    
    // dense qp
	int qp_size = d_memsize_dense_qp(&qp_dim);
	void *qp_mem = mxCalloc(qp_size,1);

	struct d_dense_qp qp;
	d_create_dense_qp(&qp_dim, &qp, qp_mem);
	d_cvt_colmaj_to_dense_qp(Hc, gc, NULL, NULL, idxb, lb, ub, Cc, lg, ug, NULL, NULL, NULL, NULL, NULL, &qp);
    
    // dense qp sol
    int qp_sol_size = d_memsize_dense_qp_sol(&qp_dim);
	void *qp_sol_mem = mxCalloc(qp_sol_size,1);

	struct d_dense_qp_sol qp_sol;
	d_create_dense_qp_sol(&qp_dim, &qp_sol, qp_sol_mem);

	// ipm arg
	int ipm_arg_size = d_memsize_dense_qp_ipm_arg(&qp_dim);
	void *ipm_arg_mem = mxCalloc(ipm_arg_size,1);

	struct d_dense_qp_ipm_arg arg;
	d_create_dense_qp_ipm_arg(&qp_dim, &arg, ipm_arg_mem);
	d_set_default_dense_qp_ipm_arg(&arg);

// 	arg.alpha_min = 1e-8;
	arg.res_g_max = 1e-4;
	arg.res_b_max = 1e-6;
	arg.res_d_max = 1e-6;
	arg.res_m_max = 1e-6;
	arg.mu0 = mu0;
	arg.iter_max = max_qp_it;
	arg.stat_max = 100;
	arg.pred_corr = pred_corr;
	arg.cond_pred_corr = cond_pred_corr;

    double *lam = (double *) mxCalloc(2*(nb+ng),sizeof(double));
    double *lam_lb = lam;
    double *lam_ub = lam+nb;
    double *lam_lg = lam+2*nb;
    double *lam_ug = lam+2*nb+ng;

	// ipm
	int ipm_size = d_memsize_dense_qp_ipm(&qp_dim, &arg);
	void *ipm_mem = mxCalloc(ipm_size,1);

	struct d_dense_qp_ipm_workspace workspace;
	d_create_dense_qp_ipm(&qp_dim, &arg, &workspace, ipm_mem);

	// call solver
	int hpipm_return = d_solve_dense_qp_ipm(&qp, &qp_sol, &arg, &workspace);

	// convert back solution
    d_cvt_dense_qp_sol_to_colmaj(&qp_sol, u, NULL, NULL, NULL, lam_lb, lam_ub, lam_lg, lam_ug, NULL, NULL);
    
    // extrac multipliers
    for(i=0;i<nb;i++)
        mu_u[i] = lam_ub[i]-lam_lb[i];
    for(i=0;i<ng;i++)
        mu_vec[i] = lam_ug[i]-lam_lg[i];
    
    //print stats
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
        
    // Free memory
    mxFree(idxb);
    mxFree(lam);
	mxFree(dense_qp_dim_mem);
	mxFree(qp_mem);
	mxFree(qp_sol_mem);
	mxFree(ipm_arg_mem);
	mxFree(ipm_mem);

	}

