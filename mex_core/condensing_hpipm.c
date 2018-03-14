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

// #include <blasfeo_target.h>
// #include <blasfeo_common.h>
// #include <blasfeo_v_aux_ext_dep.h>
// #include <blasfeo_d_aux_ext_dep.h>
// #include <blasfeo_i_aux_ext_dep.h>
// #include <blasfeo_d_aux.h>
// #include <blasfeo_d_blas.h>

#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"

#include "hpipm_d_dense_qp_dim.h"
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
// #include "hpipm_d_dense_qp_res.h"
#include "hpipm_d_cond.h"

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
    
    double *Hc = mxGetPr( mxGetField(prhs[0], 0, "Hc") );
    double *gc = mxGetPr( mxGetField(prhs[0], 0, "gc") );
    double *Cc = mxGetPr( mxGetField(prhs[0], 0, "Cc") );
    double *lbb = mxGetPr( mxGetField(prhs[0], 0, "lb_du") );
    double *ubb = mxGetPr( mxGetField(prhs[0], 0, "ub_du") );
    double *lgg = mxGetPr( mxGetField(prhs[0], 0, "lcc") );
    double *ugg = mxGetPr( mxGetField(prhs[0], 0, "ucc") );
    
    double *x = mxGetPr( mxGetField(prhs[0], 0, "dx") );
    double *u = mxGetPr( mxGetField(prhs[0], 0, "du") );
    double *pi = mxGetPr( mxGetField(prhs[0], 0, "lambda_new") );
    
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int ng = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ngN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") ); 
               
    int nb = nu;
    memcpy(x, ds0, nx*sizeof(double));
        	
	// 
	int ii, jj;

	int nb0 = nu;
    int nbN = 0;

	int nx_v[N+1];
	nx_v[0] = 0; // x0 is eliminated
	for(ii=1; ii<=N; ii++)
		nx_v[ii] = nx;

	int nu_v[N+1];
	for(ii=0; ii<N; ii++)
		nu_v[ii] = nu;
	nu_v[N] = 0;

	int nb_v[N+1];
    nb_v[0] = nb;
	for(ii=1; ii<N; ii++)
		nb_v[ii] = nb;
    nb_v[N] = 0;

	int nbu_v[N+1];
	for(ii=0; ii<N; ii++)
        nbu_v[ii] = nb;
	nbu_v[N] = 0;

	int nbx_v[N+1];
	for(ii=0; ii<=N; ii++)
        nbx_v[ii] = 0;

	int ng_v[N+1];
	for(ii=0; ii<N; ii++)
		ng_v[ii] = ng;
	ng_v[N] = ngN;

	int ns_v[N+1];
	for(ii=0; ii<=N; ii++)
		ns_v[ii] = 0;


	double b0[nx];
	dgemv_n_3l(nx, nx, 1.0, A, nx, x, 1.0, b, b0); // b0 = 1.0*b+1.0*A*x
    
	double r0[nu];
	dgemv_n_3l(nu, nx, 1.0, S, nu, x, 1.0, r, r0); // r0 = 1.0*r+1.0*S'*x

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

    // input 
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
    int *hidxb[N+1];
	double *hC[N+1];
	double *hD[N];
	double *hlg[N+1];
	double *hug[N+1];
    
    // output
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
    
    int *ptr_idx = (int *) mxCalloc((N+1)*nb,sizeof(int));
	for(ii=0; ii<=N; ii++)
		{
		hidxb[ii] = ptr_idx+ii*nb;
		for(jj=0; jj<nb_v[ii]; jj++)
			hidxb[ii][jj] = jj;
		}
	
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
    
	// ocp qp dim
	int ocp_qp_dim_size = d_memsize_ocp_qp_dim(N);
	void *ocp_qp_dim_mem = mxCalloc(ocp_qp_dim_size,1);

	struct d_ocp_qp_dim ocp_qp_dim;
	d_create_ocp_qp_dim(N, &ocp_qp_dim, ocp_qp_dim_mem);
	d_cvt_int_to_ocp_qp_dim(N, nx_v, nu_v, nbx_v, nbu_v, ng_v, ns_v, &ocp_qp_dim);
        
    // ocp qp
	int ocp_qp_size = d_memsize_ocp_qp(&ocp_qp_dim);
	void *ocp_qp_mem = mxCalloc(ocp_qp_size,1);
	struct d_ocp_qp ocp_qp;
	d_create_ocp_qp(&ocp_qp_dim, &ocp_qp, ocp_qp_mem);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, &ocp_qp);
    
    // dense qp dim
    int dense_qp_dim_size = d_memsize_dense_qp_dim();
	void *dense_qp_dim_mem = mxCalloc(dense_qp_dim_size,1);

	struct d_dense_qp_dim dense_qp_dim;
	d_create_dense_qp_dim(&dense_qp_dim, dense_qp_dim_mem);

	d_compute_qp_dim_ocp2dense(&ocp_qp_dim, &dense_qp_dim);
	
//     mexPrintf("nv=%d  ne=%d  nb=%d  ng=%d  ns=%d\n", dense_qp_dim.nv, dense_qp_dim.ne, dense_qp_dim.nb, dense_qp_dim.ng, dense_qp_dim.ns);

    // dense qp
    int dense_qp_size = d_memsize_dense_qp(&dense_qp_dim);
	void *dense_qp_mem = mxCalloc(dense_qp_size,1);
	struct d_dense_qp dense_qp;
	d_create_dense_qp(&dense_qp_dim, &dense_qp, dense_qp_mem);	
    
    // call condensing
    
    int cond_size = d_memsize_cond_qp_ocp2dense(&ocp_qp_dim);
	void *cond_mem = mxCalloc(cond_size,1);

	struct d_cond_qp_ocp2dense_workspace cond_ws;
	d_create_cond_qp_ocp2dense(&ocp_qp_dim, &cond_ws, cond_mem);
    
    // Todo: solve segfault
    d_cond_qp_ocp2dense(&ocp_qp, &dense_qp, &cond_ws); 
    d_cond_rhs_qp_ocp2dense(&ocp_qp, &dense_qp, &cond_ws);
    
    // convert
    int *idxb = (int *) mxCalloc(dense_qp_dim.nb,sizeof(int));
    d_cvt_dense_qp_to_colmaj(&dense_qp, Hc, gc, NULL, NULL, idxb, lbb, ubb, Cc, lgg, ugg, NULL, NULL, NULL, NULL, NULL);
    
    
    // Free memory

	mxFree(ptr_idx);
    mxFree(lam);
	mxFree(ocp_qp_dim_mem);
	mxFree(dense_qp_dim_mem);
    mxFree(ocp_qp_mem);
    mxFree(dense_qp_mem);
    mxFree(cond_mem);
    mxFree(idxb);
	}

