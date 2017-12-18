# MATMPC
MATMPC: MATLAB based nonlinear MPC tool

This tool aims at providing an easy-to-use nonlinear MPC implementation. The optimal control problem (OCP) that should be solved is transcribed by multiple shooting and the resulting nonlinear program (NLP) is solved by Sequential Quadratic Programming (SQP) method. Both Real-Time Iteration (RTI) and converging SQP are implemented. 

The tool supports fixed step (explicit/implicit) Runge-Kutta (RK) integrator for multiple shooting. The derivatives that are needed to perform optimization are obtained by CasADi (https://github.com/casadi/casadi/wiki), the state-of-the-art automatic/algorithmic differentiation toolbox. The Quadratic Programming (QP) problems can be solved by both dense and sparse solvers, i.e. qpOASES (https://projects.coin-or.org/qpOASES/wiki/QpoasesInstallation) and hpipm (https://github.com/giaf/hpipm). 

The most unique feature of MATMPC is that it does not require to install any external libraries. All the algorithmic routines are written in MATLAB C API and can be compiled into MEX functions using compilers that support C99 standard. In addition, CasADi and qpOASES are distributed along with MEX interface.

To use MATMPC, follow the steps below.

1. Download and install CasADi, qpOASES and hpipm (optional) by following instructions given by their developers. Note that you only need to install their MATLAB interfaces, not the entire software (e.g. unzip, run "install.m", set and save path, done!).

2. Write your own model following the styles given by examples, e.g. Inverted Pendulum, Chain of Masses.

3. In Code_generation.m, set your own sampling time and multiple shooting time.

4. Run Code_generation.m

5. In Initialization.m, set your own initialization data, e.g. initial states, referecnes and etc.

6. In Simulation.m, choose your integrator, prediction horizon and the solver. You may also need to modify the reference for constant or time-varying reference tracking problems.

7. In Draw.m, write your own plot fuctions to display your results.

8. Run Simulation.m and see the results!
