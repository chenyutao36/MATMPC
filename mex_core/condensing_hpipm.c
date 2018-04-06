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

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"

#include "hpipm_d_dense_qp_dim.h"
#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
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
    double *lgg = mxGetPr( mxGetField(prhs[0], 0, "lcc") );
    double *ugg = mxGetPr( mxGetField(prhs[0], 0, "ucc") );
    
    double *x = mxGetPr( mxGetField(prhs[0], 0, "dx") );
        
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int ng = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ngN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") ); 
                       	
	int ii, jj;

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

    // always nu
	int nb_v[N+1];
	for(ii=0; ii<N; ii++)
		nb_v[ii] = nu;
    nb_v[N] = 0;
    
    // always nu
	int nbu_v[N+1];
	for(ii=0; ii<N; ii++)
        nbu_v[ii] = nu;
	nbu_v[N] = 0;

    // always zero
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

	int *hidxb[N+1];
	for(ii=0; ii<=N; ii++)
		{
		hidxb[ii] = (int *)mxCalloc(nb_v[ii],sizeof(int));
		for(jj=0; jj<nb_v[ii]; jj++)
            hidxb[ii][jj] = jj;
		}
    
    double b0[nx];
	dgemv_n_3l(nx, nx, 1.0, A, nx, x, 1.0, b, b0);

	double r0[nu];
	dgemv_n_3l(nu, nx, 1.0, S, nu, x, 1.0, r, r0);
    
    double lg0[ng];
	dgemv_n_3l(ng, nx, -1.0, C, ng, x, 1.0, lg, lg0);

	double ug0[ng];
	dgemv_n_3l(ng, nx, -1.0, C, ng, x, 1.0, ug, ug0);

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

	for(ii=0; ii<N; ii++)
		hlb[ii] = lb+ii*nu;

	for(ii=0; ii<N; ii++)
		hub[ii] = ub+ii*nu;

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
	
    // dense qp
    int dense_qp_size = d_memsize_dense_qp(&dense_qp_dim);
	void *dense_qp_mem = mxCalloc(dense_qp_size,1);
	struct d_dense_qp dense_qp;
	d_create_dense_qp(&dense_qp_dim, &dense_qp, dense_qp_mem);	
    
    // condensing arg
    int cond_arg_size = d_memsize_cond_qp_ocp2dense_arg();
	void *cond_arg_mem = mxCalloc(cond_arg_size,1);
	struct d_cond_qp_ocp2dense_arg cond_arg;
	d_create_cond_qp_ocp2dense_arg(&cond_arg, cond_arg_mem);
	d_set_default_cond_qp_ocp2dense_arg(&cond_arg);
    
    // condensing workspace   
    int cond_size = d_memsize_cond_qp_ocp2dense(&ocp_qp_dim, &cond_arg);
	void *cond_mem = mxCalloc(cond_size,1);
	struct d_cond_qp_ocp2dense_workspace cond_ws;
	d_create_cond_qp_ocp2dense(&ocp_qp_dim, &cond_arg, &cond_ws, cond_mem);
        
    // call condensing
    d_cond_qp_ocp2dense(&ocp_qp, &dense_qp, &cond_arg, &cond_ws); 
        
    // fill in the upper triangular of H in dense_qp
    blasfeo_dtrtr_l(dense_qp_dim.nv, dense_qp.Hv, 0, 0, dense_qp.Hv, 0, 0);
        
//     int nvc = dense_qp_dim.nv;
// 	int nec = dense_qp_dim.ne;
// 	int nbc = dense_qp_dim.nb;
// 	int ngc = dense_qp_dim.ng;
// 	int nsc = dense_qp_dim.ns;

// 	printf("\nnv = %d, ne = %d, nb = %d, ng = %d, ns = %d\n\n", nvc, nec, nbc, ngc, nsc);
    
    // convert
    int idxb[dense_qp_dim.nb];
    double lbb[dense_qp_dim.nb];
    double ubb[dense_qp_dim.nb];
    
    d_cvt_dense_qp_to_colmaj(&dense_qp, Hc, gc, NULL, NULL, idxb, lbb, ubb, Cc, lgg, ugg, NULL, NULL, NULL, NULL, NULL);
    
    
    // Free memory

	for(ii=0;ii<=N;ii++){
        mxFree(hidxb[ii]);
    }
    
	mxFree(ocp_qp_dim_mem);
	mxFree(dense_qp_dim_mem);
    mxFree(ocp_qp_mem);
    mxFree(dense_qp_mem);
    mxFree(cond_arg_mem);
    mxFree(cond_mem);
	}

