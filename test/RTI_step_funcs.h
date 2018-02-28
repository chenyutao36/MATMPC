#ifndef RTI_STEP_FUNCS_H_
#define RTI_STEP_FUNCS_H_

void qp_generation(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc, double *lb_du, double *ub_du, 
        double *x0, double *x, double *u, double *y, double *yN, double *od, double *W, double *WN, double *lb, double *ub,
        double *lbN, double *ubN, double *lbu, double *ubu, int lin_obj,
        rti_step_dims *dim, rti_step_workspace *work);

void condensing(double *Q, double *S, double *R, double *A, double *B, double *Cx, double *Cu, double *CN, 
        double *gx, double *gu, double *a, double *ds0, double *lc, double *uc,
        double *G, double *Hc, double *gc, double *Cc, double *lcc, double *ucc,
        int iter, bool cond_save,
        rti_step_dims *dim, rti_step_workspace *work);

void recover(double *A, double *B, double *a, double *ds0, double *dx, double *du, rti_step_dims *dim);

void line_search(double *dx, double *du, double *x, double *u, rti_step_dims *dim);

#endif