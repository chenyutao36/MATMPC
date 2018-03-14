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

extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

// z = beta*y + alpha*A*x
void dgemv_n_3l(int m, int n, double alpha, double *A, int lda, double *x, double beta, double *y, double *z)
	{

	int ii, jj;

	double tmp;

	for(ii=0; ii<m; ii++)
		z[ii] = beta * y[ii];

	for(jj=0; jj<n; jj++)
		{
		tmp = alpha * x[jj];
		for(ii=0; ii<m; ii++)
			{
			z[ii] += A[ii+lda*jj] * tmp;
			}
		}
	
	}
 
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{
		
	// get data 
    
    double *Q = mxGetPr( mxGetField(prhs[0], 0, "Q") );
    double *S = mxGetPr( mxGetField(prhs[0], 0, "S") );
    double *R = mxGetPr( mxGetField(prhs[0], 0, "R") );
    double *A = mxGetPr( mxGetField(prhs[0], 0, "A") );
    double *B = mxGetPr( mxGetField(prhs[0], 0, "B") );
    double *C = mxGetPr( mxGetField(prhs[0], 0, "Cx") );
    double *D = mxGetPr( mxGetField(prhs[0], 0, "Cu") );
    double *CN = mxGetPr( mxGetField(prhs[0], 0, "CN") );
    double *q = mxGetPr( mxGetField(prhs[0], 0, "gx") );
    double *r = mxGetPr( mxGetField(prhs[0], 0, "gu") );   
    double *b = mxGetPr( mxGetField(prhs[0], 0, "a") );
    double *ds0 = mxGetPr( mxGetField(prhs[0], 0, "ds0") );
    double *lg = mxGetPr( mxGetField(prhs[0], 0, "lc") );
    double *ug = mxGetPr( mxGetField(prhs[0], 0, "uc") );
    double *lb = mxGetPr( mxGetField(prhs[0], 0, "lb_du") );
    double *ub = mxGetPr( mxGetField(prhs[0], 0, "ub_du") );
    
    double *x = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *pi = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    double *mu_u = mxGetPr( mxGetField(prhs[0], 0, "mu_u_new") );
    double *mu = mxGetPr( mxGetField(prhs[0], 0, "mu_new") );
    double *muN = mxGetPr( mxGetField(prhs[0], 0, "muN_new") );
    
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int ng = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ngN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") ); 
    int nbu = mxGetScalar( mxGetField(prhs[1], 0, "nbu") ); 
    double *nbu_idx = mxGetPr( mxGetField(prhs[1], 0, "nbu_idx") ); 
    
    double mu0 = mxGetScalar( mxGetField(prhs[0], 0, "mu0") ); 
    int max_qp_it = mxGetScalar( mxGetField(prhs[0], 0, "max_qp_it") );
    int pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "pred_corr") );
    int cond_pred_corr = mxGetScalar( mxGetField(prhs[0], 0, "cond_pred_corr") );
       
    int nb = nbu;
    memcpy(x, ds0, nx*sizeof(double));
        	
	// 
	int ii, jj;
// 	int i_tmp;

// 	int nb0 = nb<nu ? nb : nu;
// 	int nbN = nb-nu>0 ? nb-nu : 0;
    int nb0 = nbu;
    int nbN = 0;

    // number of states for each stage
	int nx_v[N+1];
	nx_v[0] = 0;
	for(ii=1; ii<=N; ii++)
		nx_v[ii] = nx;

    // number of controls for each stage
	int nu_v[N+1];
	for(ii=0; ii<N; ii++)
		nu_v[ii] = nu;
	nu_v[N] = 0;

    // always nbu
	int nb_v[N+1];
// 	nb_v[0] = nb<nu ? nb : nu;
    nb_v[0] = nbu;
	for(ii=1; ii<N; ii++)
// 		nb_v[ii] = nb;
        nb_v[ii] = nbu;
// 	i_tmp = nb-nu;
// 	nb_v[N] = i_tmp<0 ? 0 : i_tmp;
    nb_v[N] = 0;
    
    // always nbu
	int nbu_v[N+1];
	for(ii=0; ii<N; ii++)
// 		nbu_v[ii] = nb<nu ? nb : nu;
        nbu_v[ii] = nbu;
	nbu_v[N] = 0;

    // always zero
	int nbx_v[N+1];
	for(ii=0; ii<=N; ii++)
		{
// 		i_tmp = nb_v[ii]-nbu_v[ii];
// 		nbx_v[ii] = i_tmp>=0 ? i_tmp : 0;
        nbx_v[ii] = 0;
		}

	int ng_v[N+1];
	for(ii=0; ii<N; ii++)
		ng_v[ii] = ng;
	ng_v[N] = ngN;

	int ns_v[N+1];
	for(ii=0; ii<=N; ii++)
		ns_v[ii] = 0;

	int *hidxb[N+1];
	int *ptr_idx = (int *) mxCalloc((N+1)*nb, sizeof(int));
	for(ii=0; ii<=N; ii++)
		{
		hidxb[ii] = ptr_idx+ii*nb;
		for(jj=0; jj<nb_v[ii]; jj++)
			hidxb[ii][jj] = (int) nbu_idx[jj]-1;
		}

	double b0[nx];
	dgemv_n_3l(nx, nx, 1.0, A, nx, x, 1.0, b, b0); // b0 = 1.0*b+1.0*A*x
    
	double r0[nu];
	dgemv_n_3l(nu, nx, 1.0, S, nu, x, 1.0, r, r0); // r0 = 1.0*r+1.0*S*x

	double lb0[nb0];
	for(ii=0; ii<nb0; ii++)
		lb0[ii] = lb[ii];

	double ub0[nb0];
	for(ii=0; ii<nb0; ii++)
		ub0[ii] = ub[ii];

	double lbN[nbN];
	double ubN[nbN];

	double lg0[ng];
	dgemv_n_3l(ng, nx, -1.0, C, ng, x, 1.0, lg, lg0); // lg0 = 1.0*lg-1.0*C*x

	double ug0[ng]; 
	dgemv_n_3l(ng, nx, -1.0, C, ng, x, 1.0, ug, ug0); // ug0 = 1.0*ug-1.0*C*x


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

    for(ii=1; ii<N; ii++)
		hA[ii] = A+ii*nx*nx;

	for(ii=0; ii<N; ii++)
		hB[ii] = B+ii*nx*nu;

	hb[0] = b0;
	for(ii=1; ii<N; ii++)
		hb[ii] = b+ii*nx;

	for(ii=1; ii<=N; ii++)
		hQ[ii] = Q+ii*nx*nx;
        
	for(ii=1; ii<N; ii++)
		hS[ii] = S+ii*nx*nu;

    for(ii=0; ii<N; ii++)
		hR[ii] = R+ii*nu*nu;

	for(ii=1; ii<=N; ii++)
		hq[ii] = q+ii*nx;
		
	hr[0] = r0;
	for(ii=1; ii<N; ii++)
		hr[ii] = r+ii*nu;

	for(ii=0; ii<nbN; ii++)
		lbN[ii] = lb[nu+ii];
	hlb[0] = lb0;
	for(ii=1; ii<N; ii++)
		hlb[ii] = lb+ii*nb;
	hlb[N] = lbN;

	for(ii=0; ii<nbN; ii++)
		ubN[ii] = ub[nu+ii];
	hub[0] = ub0;
	for(ii=1; ii<N; ii++)
		hub[ii] = ub+ii*nb;
	hub[N] = ubN;

	for(ii=1; ii<N; ii++)
		hC[ii] = C+ii*ng*nx;
	hC[N] = CN;

	for(ii=0; ii<N; ii++)
		hD[ii] = D+ii*ng*nu;

	hlg[0] = lg0;
	for(ii=1; ii<=N; ii++)
		hlg[ii] = lg+ii*ng;

	hug[0] = ug0;
	for(ii=1; ii<=N; ii++)
		hug[ii] = ug+ii*ng;
	
	for(ii=0; ii<=N; ii++)
		hx[ii] = x+ii*nx;

	for(ii=0; ii<N; ii++)
		hu[ii] = u+ii*nu;
	
	for(ii=0; ii<N; ii++)
		hpi[ii] = pi+(ii+1)*nx;
    
    double *lam = (double *) mxCalloc((2*(N+1)*nb+2*N*ng+2*ngN),sizeof(double));
    for(ii=0; ii<=N; ii++)
        {
		hlam_lb[ii] = lam+ii*nb;
		hlam_ub[ii] = lam+(N+1)*nb+ii*nb;
		hlam_lg[ii] = lam+2*(N+1)*nb+ii*ng;
        hlam_ug[ii] = lam+2*(N+1)*nb+N*ng+ngN+ii*ng;
        }
    
	// qp dim
	int dim_size = d_memsize_ocp_qp_dim(N);
// 	void *dim_mem = mxMalloc(dim_size);
    void *dim_mem = mxCalloc(dim_size,1);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);
	d_cvt_int_to_ocp_qp_dim(N, nx_v, nu_v, nbx_v, nbu_v, ng_v, ns_v, &dim);

	// qp
	int qp_size = d_memsize_ocp_qp(&dim);
// 	void *qp_mem = mxMalloc(qp_size);
    void *qp_mem = mxCalloc(qp_size,1);
	struct d_ocp_qp qp;
	d_create_ocp_qp(&dim, &qp, qp_mem);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, &qp);


	// qp sol
	int qp_sol_size = d_memsize_ocp_qp_sol(&dim);
// 	void *qp_sol_mem = mxMalloc(qp_sol_size);
    void *qp_sol_mem = mxCalloc(qp_sol_size,1);
	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(&dim, &qp_sol, qp_sol_mem);


	// ipm arg
	int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(&dim);
// 	void *ipm_arg_mem = mxMalloc(ipm_arg_size);
    void *ipm_arg_mem = mxCalloc(ipm_arg_size,1);

	struct d_ocp_qp_ipm_arg arg;
	d_create_ocp_qp_ipm_arg(&dim, &arg, ipm_arg_mem);
	d_set_default_ocp_qp_ipm_arg(&arg);

	arg.alpha_min = 1e-10;
	arg.res_g_max = 1e-4;
	arg.res_b_max = 1e-6;
	arg.res_d_max = 1e-6;
	arg.res_m_max = 1e-6;
	arg.mu0 = mu0;
	arg.iter_max = max_qp_it;
	arg.stat_max = 100;
	arg.pred_corr = pred_corr;
	arg.cond_pred_corr = cond_pred_corr;


	// ipm
	int ipm_size = d_memsize_ocp_qp_ipm(&dim, &arg);
// 	void *ipm_mem = mxMalloc(ipm_size);
    void *ipm_mem = mxCalloc(ipm_size,1);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim, &arg, &workspace, ipm_mem);

	// call solver
	int hpipm_return = d_solve_ocp_qp_ipm(&qp, &qp_sol, &arg, &workspace);

	// convert back solution
	d_cvt_ocp_qp_sol_to_colmaj(&qp_sol, hu, hx, NULL, NULL, hpi, hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

    
    // extrac multipliers
    for(ii=0;ii<N;ii++){
        for(jj=0;jj<nb;jj++)
            mu_u[ii*nb+jj] = hlam_ub[ii][jj] - hlam_lb[ii][jj];
        
        for(jj=0;jj<ng;jj++)
            mu[ii*ng+jj] = hlam_ug[ii][jj] - hlam_lg[ii][jj];
    }
    for(jj=0;jj<ngN;jj++)
        muN[jj] = hlam_ug[N][jj] - hlam_lg[N][jj]; 
            
    char *nTrans = "N", *Trans="T";
    double one_d = 1.0, minus_one_d = -1.0;
    int one_i = 1;
    memcpy(pi, q, nx*sizeof(double));
    dgemv_(Trans,&nx,&nx,&minus_one_d,A,&nx,pi+nx,&one_i,&minus_one_d,pi,&one_i);
    dgemv_(nTrans,&nx,&nx,&minus_one_d,Q,&nx,x,&one_i,&one_d,pi,&one_i);
    dgemv_(nTrans,&nx,&nu,&minus_one_d,S,&nx,u,&one_i,&one_d,pi,&one_i);
    dgemv_(Trans,&ng,&nx,&minus_one_d,C,&ng,mu,&one_i,&one_d,pi,&one_i);
    
    // print stats
//     mexPrintf("res_g:%5.3e  res_b:%5.3e  res_d:%5.3e  res_m:%5.3e", workspace.qp_res[0],workspace.qp_res[1],workspace.qp_res[2],workspace.qp_res[3]);
//     mexPrintf("  No. of It: %d\n", workspace.iter);
    
    // Free memory
	mxFree(ptr_idx);
    mxFree(lam);
	mxFree(dim_mem);
	mxFree(qp_mem);
	mxFree(qp_sol_mem);
	mxFree(ipm_arg_mem);
	mxFree(ipm_mem);

	}

