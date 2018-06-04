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
#include "hpipm_d_part_cond.h"
 
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
    
    int nx = mxGetScalar( mxGetField(prhs[1], 0, "nx") );
    int nu = mxGetScalar( mxGetField(prhs[1], 0, "nu") );
    int ng = mxGetScalar( mxGetField(prhs[1], 0, "nc") ); 
    int ngN = mxGetScalar( mxGetField(prhs[1], 0, "ncN") );
    int N = mxGetScalar( mxGetField(prhs[1], 0, "N") );
    int N2 = mxGetScalar( mxGetField(prhs[1], 0, "N2") ); 
    int nbx = mxGetScalar( mxGetField(prhs[1], 0, "nbx") );
    double *nbx_idx = mxGetPr( mxGetField(prhs[1], 0, "nbx_idx") );
                   	 
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
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &ocp_qp);

    // partial ocp qp dim
    int nx2[N2+1];
	int nu2[N2+1];
	int nb2[N2+1];
	int nbx2[N2+1];
	int nbu2[N2+1];
	int ng2[N2+1];
	int ns2[N2+1];
        
    int ocp_qp_dim_size2 = d_memsize_ocp_qp_dim(N2);
	void *ocp_qp_dim_mem2 = mxCalloc(ocp_qp_dim_size2,1);

	struct d_ocp_qp_dim ocp_qp_dim2;
	d_create_ocp_qp_dim(N2, &ocp_qp_dim2, ocp_qp_dim_mem2);
	d_cvt_int_to_ocp_qp_dim(N2, nx2, nu2, nbx2, nbu2, ng2, ns2, &ocp_qp_dim2);
    
    int block_size[N2+1];
    
    d_compute_block_size_cond_qp_ocp2ocp(N, N2, block_size);
    d_compute_qp_dim_ocp2ocp(&ocp_qp_dim, block_size, &ocp_qp_dim2);
        
    // partial ocp qp
    int ocp_qp_size2 = d_memsize_ocp_qp(&ocp_qp_dim2);
	void *ocp_qp_mem2 = mxCalloc(ocp_qp_size2,1);
	struct d_ocp_qp ocp_qp2;
	d_create_ocp_qp(&ocp_qp_dim2, &ocp_qp2, ocp_qp_mem2);
    
    // partial condensing arg
    int part_cond_arg_size = d_memsize_cond_qp_ocp2ocp_arg(N2);
    void *part_cond_arg_mem = mxCalloc(part_cond_arg_size,1);
    struct d_cond_qp_ocp2ocp_arg part_cond_arg;
    d_create_cond_qp_ocp2ocp_arg(N2, &part_cond_arg, part_cond_arg_mem);
    d_set_default_cond_qp_ocp2ocp_arg(N2, &part_cond_arg);

    // partial condensing workspace
	int part_cond_size = d_memsize_cond_qp_ocp2ocp(&ocp_qp_dim, block_size, &ocp_qp_dim2, &part_cond_arg);
	void *part_cond_mem = mxCalloc(part_cond_size,1);
	struct d_cond_qp_ocp2ocp_workspace part_cond_ws;
	d_create_cond_qp_ocp2ocp(&ocp_qp_dim, block_size, &ocp_qp_dim2, &part_cond_arg, &part_cond_ws, part_cond_mem);
    
    // call partial condensing
    d_cond_qp_ocp2ocp(&ocp_qp, &ocp_qp2, &part_cond_arg, &part_cond_ws);
    
//     for (ii=0;ii<=N2;ii++){
//         mexPrintf("nx2[%d]=%d   nu2[%d]=%d   nb2[%d]=%d   nbx2[%d]=%d   nbu2[%d]=%d   ng2[%d]=%d\n", ii,ocp_qp_dim2.nx[ii],ii,ocp_qp_dim2.nu[ii],
//                 ii,ocp_qp_dim2.nb[ii],ii,ocp_qp_dim2.nbx[ii],ii,ocp_qp_dim2.nbu[ii],ii,ocp_qp_dim2.ng[ii]);
//     }
    
	// convert partial condensed QP to colmaj
    plhs[0] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], ocp_qp_dim2.nx[0]*(N2+1), mxREAL);  
    double *Qp = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(ocp_qp_dim2.nu[0], ocp_qp_dim2.nx[0]*N2, mxREAL);
    double *Sp = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(ocp_qp_dim2.nu[0], ocp_qp_dim2.nu[0]*N2, mxREAL);  
    double *Rp = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], N2+1, mxREAL);  
    double *qp = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(ocp_qp_dim2.nu[0], N2, mxREAL);  
    double *rp = mxGetPr(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], ocp_qp_dim2.nx[0]*N2, mxREAL);  
    double *Ap = mxGetPr(plhs[5]);
    plhs[6] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], ocp_qp_dim2.nu[0]*N2, mxREAL);  
    double *Bp = mxGetPr(plhs[6]);
    plhs[7] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], N2, mxREAL);  
    double *bp = mxGetPr(plhs[7]);
    plhs[8] = mxCreateDoubleMatrix(ocp_qp_dim2.nbx[1]*N2, 1, mxREAL);  
    double *ubxp = mxGetPr(plhs[8]);
    plhs[9] = mxCreateDoubleMatrix(ocp_qp_dim2.nbx[1]*N2, 1, mxREAL);  
    double *lbxp = mxGetPr(plhs[9]);
    plhs[10] = mxCreateDoubleMatrix(ocp_qp_dim2.nbu[0]*N2, 1, mxREAL);  
    double *ubup = mxGetPr(plhs[10]);
    plhs[11] = mxCreateDoubleMatrix(ocp_qp_dim2.nbu[0]*N2, 1, mxREAL);  
    double *lbup = mxGetPr(plhs[11]);
    
    plhs[12] = mxCreateDoubleMatrix(ocp_qp_dim2.ng[0], ocp_qp_dim2.nx[0]*N2, mxREAL);  
    double *Cp = mxGetPr(plhs[12]);
    plhs[13] = mxCreateDoubleMatrix(ocp_qp_dim2.ng[0], ocp_qp_dim2.nu[0]*N2, mxREAL);  
    double *Dp = mxGetPr(plhs[13]);
    plhs[14] = mxCreateDoubleMatrix(ocp_qp_dim2.ng[N2], ocp_qp_dim2.nx[N2], mxREAL);  
    double *CNp = mxGetPr(plhs[14]);
    plhs[15] = mxCreateDoubleMatrix(ocp_qp_dim2.ng[0]*N2+ocp_qp_dim2.ng[N2], 1, mxREAL);  
    double *ugp = mxGetPr(plhs[15]);
    plhs[16] = mxCreateDoubleMatrix(ocp_qp_dim2.ng[0]*N2+ocp_qp_dim2.ng[N2], 1, mxREAL);  
    double *lgp = mxGetPr(plhs[16]);
    
    plhs[17] = mxCreateDoubleMatrix(ocp_qp_dim2.nx[0], 1, mxREAL);  
    double *ds0p = mxGetPr(plhs[17]);
    d_cvt_ocp_qp_to_colmaj_ubx(0, &ocp_qp2, ds0p);
    
//     plhs[18] = mxCreateDoubleScalar(ocp_qp_dim2.nx[0]); 
//     plhs[19] = mxCreateDoubleScalar(ocp_qp_dim2.nu[0]); 
//     plhs[20] = mxCreateDoubleScalar(ocp_qp_dim2.nbx[1]); 
//     plhs[21] = mxCreateDoubleScalar(ocp_qp_dim2.nbu[0]); 
//     plhs[22] = mxCreateDoubleScalar(ocp_qp_dim2.ng[0]); 
//     plhs[23] = mxCreateDoubleScalar(ocp_qp_dim2.ng[N2]); 
        
    for (ii=0;ii<=N2;ii++){
        d_cvt_ocp_qp_to_colmaj_Q(ii, &ocp_qp2, Qp+ii*ocp_qp_dim2.nx[0]*ocp_qp_dim2.nx[0]);
        d_cvt_ocp_qp_to_colmaj_q(ii, &ocp_qp2, qp+ii*ocp_qp_dim2.nx[0]);
        
        d_cvt_ocp_qp_to_colmaj_ug(ii, &ocp_qp2, ugp+ii*ocp_qp_dim2.ng[0]);
        d_cvt_ocp_qp_to_colmaj_lg(ii, &ocp_qp2, lgp+ii*ocp_qp_dim2.ng[0]);
    }
    for (ii=0;ii<N2;ii++){
        d_cvt_ocp_qp_to_colmaj_S(ii, &ocp_qp2, Sp+ii*ocp_qp_dim2.nu[0]*ocp_qp_dim2.nx[0]);
        d_cvt_ocp_qp_to_colmaj_R(ii, &ocp_qp2, Rp+ii*ocp_qp_dim2.nu[0]*ocp_qp_dim2.nu[0]);  
        d_cvt_ocp_qp_to_colmaj_r(ii, &ocp_qp2, rp+ii*ocp_qp_dim2.nu[0]);
        d_cvt_ocp_qp_to_colmaj_A(ii, &ocp_qp2, Ap+ii*ocp_qp_dim2.nx[0]*ocp_qp_dim2.nx[0]);  
        d_cvt_ocp_qp_to_colmaj_B(ii, &ocp_qp2, Bp+ii*ocp_qp_dim2.nx[0]*ocp_qp_dim2.nu[0]); 
        d_cvt_ocp_qp_to_colmaj_b(ii, &ocp_qp2, bp+ii*ocp_qp_dim2.nx[0]);
        d_cvt_ocp_qp_to_colmaj_ubx(ii+1, &ocp_qp2, ubxp+ii*ocp_qp_dim2.nbx[1]);
        d_cvt_ocp_qp_to_colmaj_lbx(ii+1, &ocp_qp2, lbxp+ii*ocp_qp_dim2.nbx[1]);
        d_cvt_ocp_qp_to_colmaj_ubu(ii, &ocp_qp2, ubup+ii*ocp_qp_dim2.nbu[0]);
        d_cvt_ocp_qp_to_colmaj_lbu(ii, &ocp_qp2, lbup+ii*ocp_qp_dim2.nbu[0]);
        
        d_cvt_ocp_qp_to_colmaj_C(ii, &ocp_qp2, Cp+ii*ocp_qp_dim2.ng[0]*ocp_qp_dim2.nx[0]);
        d_cvt_ocp_qp_to_colmaj_D(ii, &ocp_qp2, Dp+ii*ocp_qp_dim2.ng[0]*ocp_qp_dim2.nu[0]); 
        
    }
    d_cvt_ocp_qp_to_colmaj_C(N2, &ocp_qp2, CNp);
    
    // Free memory
	for(ii=0;ii<=N;ii++){
        mxFree(hidxb[ii]);
        mxFree(hlb[ii]);
        mxFree(hub[ii]);
    }
    
	mxFree(ocp_qp_dim_mem);
	mxFree(ocp_qp_mem);
    mxFree(ocp_qp_dim_mem2);
    mxFree(ocp_qp_mem2);
    
    mxFree(part_cond_arg_mem);
    mxFree(part_cond_mem);
 
   
}

