# MATMPC
MATMPC: MATLAB based nonlinear MPC tool

This tool aims at providing an easy-to-use nonlinear MPC implementation. The optimal control problem (OCP) that should be solved is transcribed by multiple shooting and the resulting nonlinear program (NLP) is solved by Sequential Quadratic Programming (SQP) method.

The tool supports fixed step (explicit/implicit) Runge-Kutta (RK) integrator for multiple shooting. The derivatives that are needed to perform optimization are obtained by CasADi (https://github.com/casadi/casadi/wiki), the state-of-the-art automatic/algorithmic differentiation toolbox. The Quadratic Programming (QP) problems can be solved by both dense and sparse solvers, i.e. qpOASES (https://projects.coin-or.org/qpOASES/wiki/QpoasesInstallation) and hpipm (https://github.com/giaf/hpipm). 

The most unique feature of MATMPC is that it does not require to install any external libraries. Execpt external QP solvers, all the algorithmic routines are written directly using MATLAB C API and can be compiled into MEX functions using compilers that belong to GCC class (e.g. GCC and MinGW). MATMPC employs MATLAB built-in linear algebra library originated from BLAS and LAPACK. 

To use MATMPC, follow the steps below.

1. Download and install CasADi, qpOASES and hpipm (optional) by following instructions given by their developers. MATMPC currently is compatible with CasADi 3.3.0 and qpOASES 3.2.1

2. Write your own model using following the styles given by examples, e.g. Inverted Pendulum, Chain of Masses.

3. In the model file you created, set your own sampling time and multiple shooting time. In Code_generation.m, set the number of integration steps per shooting interval.

4. Run Code_generation.m

5. In InitData.m, set your own initialization data, e.g. data path, initial states, referecnes and etc.

6. In Simulation.m, choose your integrator, set the prediction horizon and the solver options. You may also need to modify the reference type for constant or time-varying reference tracking problems.

7. In Draw.m, write your own plot fuctions to display your results.

8. Run Simulation.m and see the results!
