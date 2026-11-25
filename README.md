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

## Prediction grid and move blocking

Set `opt.nonuniform_grid = 1` in `Simulation.m` to use unequal shooting intervals. Set `settings.r` to the interval count and optionally set `settings.index_T` to the interval boundaries in units of `settings.Ts_st`. For example, `[0 1 3 6]` gives three intervals of lengths one, two, and three base sampling periods. The prediction horizon passed to `InitMemory` uses `settings.N = settings.r` in this mode. Reference preview follows the same boundaries.

For uniform shooting with move blocking, select `qpoases_mb`. Set `settings.r` to the number of control blocks and enter `r + 1` strictly increasing integer boundaries from zero to `settings.N`. The default boundaries cover 80 uniform steps in 10 blocks. Nonuniform shooting and move blocking are separate solver modes.

`mpc_nmpcsolver` reports step acceptance in `output.info.lineSearchAccepted`. When the line search rejects all trial steps, the solver returns the last accepted iterate.

## Simulink examples

Build the inverted pendulum and initialize the controller workspace before opening either supplied model.

```matlab
Build_Model('InvertedPendulum')
Initialization_Simulink
open_system('MATMPC_SIMULINK')
```

Use `open_system('MATMPC_SIMULINK_R2017b')` for the legacy model.
Run `addpath('tests'); smoke_simulink_inverted_pendulum` to build and check both models.

## Shared C kernels

The explicit Runge–Kutta integrator, default full condensing, and condensed QP solution recovery live in `core/`. The existing MEX entry points link these same sources. The C APIs operate on ordinary arrays and caller-owned memory and can be built independently of MATLAB, CasADi, and the QP solvers.

```sh
cmake -S core -B build/shared-core
cmake --build build/shared-core --config Release
cmake -E chdir build/shared-core ctest -C Release --output-on-failure
```

Use `-DBUILD_TESTING=OFF` when building the library alone. The ERK caller supplies model callbacks and the queried scratch storage. Recovery reads the QP data and writes into separate trajectory and multiplier buffers. The standalone build uses column-major C linear algebra for the condensing kernels. MATLAB uses an adapter to its existing BLAS. The MATLAB setup above rebuilds the MEX gateways when using these sources.

Full condensing exposes the gradient and bounds separately from the Hessian and constraint matrices in C. In MATLAB, `qp = Condensing(mem, settings)` returns all eight condensed arrays for inspection or comparison.

## MATLAB controller sessions

Use a controller object to manage one solve sequence and inspect its numeric state in the workspace.

```matlab
settings = Build_Model('InvertedPendulum');
controller = MATMPCController(settings);
initial = controller.snapshot();
[control, solution] = controller.step(initial.input.x0);
debug = controller.snapshot(true);
controller.reset();
delete(controller);
```

`step` accepts a second argument with updated references, weights, bounds, or online model parameters. `snapshot()` returns ordinary MATLAB data without copying QP work arrays. `snapshot(true)` also includes the latest available QP arrays and can be saved with `save`. Controllers own their backend warm starts. Multiple live controllers use the same active generated model. Release the instances before rebuilding the MEX files, calling `clear mex`, or starting an example that clears the MEX state. Scripts and Simulink can continue to call the function interface.

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
