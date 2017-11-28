#include "mex.h"
#include <qpOASES_e.h>
#include "assert.h"

static void *raw_memory_ptr = NULL;
static double *xOpt = NULL, *yOpt = NULL;
static int qp_instance =0;

void ExitFcn(){
    if (raw_memory_ptr!=NULL){
        mxFree(raw_memory_ptr);
        mxFree(xOpt);
        mxFree(yOpt);
    }
}

void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{ 
	USING_NAMESPACE_QPOASES

	double stat, feas, cmpl;

	/* Setup data of first QP. */
	double H[2*2] = { 1.0, 0.0, 0.0, 0.5 };
	double A[1*2] = { 1.0, 1.0 };
	double g[2] = { 1.5, 1.0 };
	double lb[2] = { 0.5, -2.0 };
	double ub[2] = { 5.0, 2.0 };
	double lbA[1] = { -1.0 };
	double ubA[1] = { 2.0 };

	/* Setup data of second QP. */
	double g_new[2] = { 1.0, 1.5 };
	double lb_new[2] = { 0.0, -1.0 };
	double ub_new[2] = { 5.0, -0.5 };
	double lbA_new[1] = { -2.0 };
	double ubA_new[1] = { 1.0 };
    


	/* Setting up QProblem object. */
	Options options;

    QProblem *example;
    
    int nWSR;
     
	double objVal;
    
    if(raw_memory_ptr==NULL){
        
        int memory_size = QProblem_calculateMemorySize(2, 1);
        raw_memory_ptr = mxMalloc(memory_size);
        char *ptr_end =  QProblem_assignMemory(2, 1, &example, raw_memory_ptr);
        assert((char*)raw_memory_ptr + memory_size >= ptr_end); (void) ptr_end;
        mexMakeMemoryPersistent(raw_memory_ptr);
        
        xOpt = (double *)mxCalloc(2,sizeof(double));
        mexMakeMemoryPersistent(xOpt);
        
        yOpt = (double *)mxCalloc(2+1,sizeof(double));
        mexMakeMemoryPersistent(yOpt);
        
        QProblemCON( example, 2, 1, HST_UNKNOWN );
        Options_setToDefault( &options );
        QProblem_setOptions( example,options );

        /* Solve first QP. */

        QProblem_setOptions( example,options );

        nWSR = 10;
        QProblem_init( example,H,g,A,lb,ub,lbA,ubA, &nWSR,0 );

        /* Get and print solution of first QP. */
        QProblem_getPrimalSolution( example,xOpt );
        QProblem_getDualSolution(   example,yOpt );
        objVal = QProblem_getObjVal( example );
        printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
                xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2], objVal );

        qpOASES_getKktViolation( 2,1, H,g,A,lb,ub,lbA,ubA, xOpt,yOpt, &stat,&feas,&cmpl );
        printf("KKT violations:\n\n");
        printf("stat = %e, feas = %e, cmpl = %e\n\n", stat, feas, cmpl);
        
        mexAtExit(ExitFcn);
    }
    
    char *ptr_end =  QProblem_assignMemory(2, 1, &example, raw_memory_ptr);	   
    
	/* Solve second QP. */
    
	nWSR = 10;
	QProblem_hotstart( example,g_new,lb_new,ub_new,lbA_new,ubA_new, &nWSR,0 );

	/* Get and print solution of second QP. */
	QProblem_getPrimalSolution( example,xOpt );
	QProblem_getDualSolution(   example,yOpt );
	objVal = QProblem_getObjVal( example );
	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",
			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2], objVal );

	qpOASES_getKktViolation( 2,1, H,g_new,A,lb_new,ub_new,lbA_new,ubA_new, xOpt,yOpt, &stat,&feas,&cmpl );
	printf("KKT violations:\n\n");
	printf("stat = %e, feas = %e, cmpl = %e\n\n", stat, feas, cmpl);
    
}