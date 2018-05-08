#ifndef CONDENSING_BLASFEO_H_
#define CONDENSING_BLASFEO_H_

typedef struct
{

    struct blasfeo_dmat *StQ;
    struct blasfeo_dmat *BtAt;
    struct blasfeo_dmat *HtWt;
    struct blasfeo_dmat *Gt;
    struct blasfeo_dmat *Rt;
    struct blasfeo_dmat *Cgx_dmat;
    struct blasfeo_dmat *CgN_dmat;
    struct blasfeo_dmat *Cx_dmat;
    struct blasfeo_dmat *Hc_dmat;
    struct blasfeo_dmat *Ccg_dmat;
    struct blasfeo_dmat *Ccx_dmat;
        
    struct blasfeo_dvec *rq;
    struct blasfeo_dvec *a_dvec;
    struct blasfeo_dvec *L;
    struct blasfeo_dvec *gcw;
    struct blasfeo_dvec *bg;
    struct blasfeo_dvec *bx;

} Condensing_Blasfeo_workspace;

int Condensing_Blasfeo_workspace_calculate_size(int nx, int nu, int N, int nc, int ncN, int nbx);

Condensing_Blasfeo_workspace* Condensing_Blasfeo_workspace_cast(int nx, int nu, int N, int nc, int ncN, int nbx, void *raw_memory);

int align_char_to(int num, char **c_ptr);

void assign_blasfeo_dmat_mem(int m, int n, struct blasfeo_dmat *sA, char **ptr);

void assign_blasfeo_dvec_mem(int n, struct blasfeo_dvec *sv, char **ptr);

#endif