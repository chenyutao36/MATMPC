#ifndef RTI_STEP_FUNCS_H_
#define RTI_STEP_FUNCS_H_

void qp_generation(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CxN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc, double *lb_du, double *ub_du, 
        double *x0, double *z, double *xN, double *y, double *yN, double *od, double *W, double *WN, double *lb, double *ub,
        double *lbN, double *ubN, double *lbu, double *ubu, int lin_obj,
        rti_step_dims *dim, rti_step_workspace *work);

void condensing(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CxN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc,
        double *G, double *Hc, double *gc, double *Cc, double *lcc, double *ucc,
        int iter, bool cond_save,
        rti_step_dims *dim, rti_step_workspace *work);

void recover(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *CxN,
        double *gx, double *a, double *ds0, double *du, double *multipliers,
        double *dz, double *dxN, double *lambda_new, double *mu_new, double *muN_new, double *mu_u_new,
        rti_step_dims *dim);

void line_search(double *dz, double *dxN, double *lambda_new, double *mu_new, double *muN_new, 
        double *z, double *xN, double *lambda, double *mu, double *muN,
        rti_step_dims *dim);

#endif