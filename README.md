# MATMPC

MATMPC solves nonlinear model predictive control problems in MATLAB. It transcribes the optimal control problem by multiple shooting and solves the resulting nonlinear program with sequential quadratic programming. It uses CasADi for derivatives, fixed-step explicit and implicit Runge–Kutta integrators for multiple shooting, and MATLAB's BLAS and LAPACK libraries for the default MEX routines. See the [MATMPC documentation](https://matmpc-docs.readthedocs.io/en/latest/) for the mathematics and tuning parameters.

## Setup

1. Install the [CasADi 3.3.0 MATLAB package](https://github.com/casadi/casadi/releases/tag/3.3.0) for your platform, extract it, and add its top-level directory to the MATLAB path.
2. Install a C compiler listed for your [MATLAB release and platform](https://www.mathworks.com/support/requirements/supported-compilers.html). Select it in MATLAB with `mex -setup C`. Both model and core MEX files use that selection.
3. Check your default solver installation. MATMPC ships qpOASES MEX binaries for Windows x64, Linux x64, and Intel macOS. A native Apple Silicon build requires a compatible CasADi package and a qpOASES MEX file matching the extension reported by `mexext`. The CasADi 3.3.0 release provides an Intel macOS package.

Use the following commands from the repository root.

```matlab
addpath('/path/to/extracted/casadi')
mex -setup C
Build_Model('InvertedPendulum')
Simulation
```

`Build_Model` generates the model and compiles its MEX files without prompting. A successful build writes `data/settings.mat`. `Model_Generation` is the interactive entry point. Set its `model` variable to select an example or another model script on the MATLAB path. Use `Build_Model('InvertedPendulum', false)` to generate C sources only. For a standalone rebuild of the core from previously generated C sources, run `addpath('mex_core'); Compile_Mex` from the repository root. Use `Build_Model` to rebuild the model MEX files and update the saved simulation settings together.

To run your own model, follow an [example](examples), configure its shooting interval and initial state, references, and bounds in `InitData.m`, choose your horizon and solver settings in `Simulation.m`, and run the same build and simulation steps. Plotting is in `Draw.m`.

## Smoke test

The test builds the shipped inverted pendulum and runs five closed-loop samples with the default qpOASES solver. It checks that the solver returns finite outputs, executes an SQP iteration each step, and keeps the applied controls and observed position within the example's bounds.

```sh
matlab -batch "addpath('/path/to/extracted/casadi'); addpath('tests'); smoke_inverted_pendulum"
```

The test builds into the repository's usual model and MEX directories. Run it in a clean checkout to keep any other generated model from being replaced.

## Optional solvers

MATMPC also has interfaces for [Ipopt](https://github.com/coin-or/Ipopt), [HPIPM](https://github.com/giaf/hpipm), [OSQP](https://osqp.org/), and MATLAB's `quadprog`. See the [HPIPM installation instructions](doc/HPIPM-Tutorial.txt) for `hpipm_sparse` and `hpipm_pcond`.

## Open in MATLAB Online

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=chenyutao36/MATMPC)
