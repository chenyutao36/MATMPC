#include "mex.h"
#include "sim.h"


sim_opts* sim_opts_create(const mxArray *mem)
{
    sim_opts* opts = (sim_opts*)mxMalloc(sizeof(sim_opts));
    mexMakeMemoryPersistent(opts);
       
    opts->nx = mxGetScalar( mxGetField(mem, 0, "nx") );
    opts->nu = mxGetScalar( mxGetField(mem, 0, "nu") );
    opts->nz = mxGetScalar( mxGetField(mem, 0, "nz") );
    opts->num_stages = mxGetScalar( mxGetField(mem, 0, "num_stages") );
    opts->num_steps = mxGetScalar( mxGetField(mem, 0, "num_steps") );
    opts->h = mxGetScalar( mxGetField(mem, 0, "h") );
    
	return opts;
}

sim_in* sim_in_create(sim_opts *opts)
{
    sim_in* in = (sim_in*)mxMalloc(sizeof(sim_in));
    mexMakeMemoryPersistent(in);
    
    return in;
}

sim_out* sim_out_create(sim_opts *opts)
{
    sim_out* out = (sim_out*)mxMalloc(sizeof(sim_out));
    mexMakeMemoryPersistent(out);
    
    return out;
}

void sim_opts_free(sim_opts *opts)
{
    mxFree(opts);
}

void sim_in_free(sim_in *in)
{
    mxFree(in);
}

void sim_out_free(sim_out *out)
{
    mxFree(out);
}