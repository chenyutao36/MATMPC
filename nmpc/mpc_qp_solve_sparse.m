function [dz, dxN, lambda, mu, muN, cpt_qp] = mpc_qp_solve_sparse(Q_h,S,R,A,B,Cx,Cu,gx,gu,a,ds0,lc,uc,lb_du,ub_du,CxN,sizes,mem)

nx = sizes.nx;
nu = sizes.nu;
N = sizes.N;

% nb  = nu+nx/2;		% number of two-sided box constraints
nb = nu;
ng  = sizes.nc;            % number of two-sided general constraints
% ngN = nx;           % number of two-sided general constraints on last stage
ngN = sizes.ncN;

time_invariant = 0; % time_invariant (0) vs time_variant (1) problems

nbu = min(nb, nu);

AA = A;
BB = B;
bb = a;

% cost function

QQ = Q_h(:,1:N*nx);
Qf = Q_h(:,N*nx:end);
SS = S;
RR = R;
qq = gx(:,1:N);
qf = gx(:,N+1);
rr = gu;

% general constraints

CN = CxN;
%if(ng>0)
%	if(nb==0)
%		for ii=1:min(nu,ng)
%			D(ii,ii) = 1.0;
%		end
%		for ii=min(nu,ng)+1:ng
%			C(ii,ii-nu) = 1.0;
%		end
%	else
%		for ii=1:ng
%			C(ii,ii) = 1.0;
%		end
%	end
%end
% if(ngN>0)
% 	for ii=1:ngN
% 		CN(ii,ii) = 1.0;
% 	end
% end
CC = Cx;
DD = Cu;

% box constraints
llb = [reshape(lb_du,[nu,N]), lb_du(1:nu)];
uub = [reshape(ub_du,[nu,N]), ub_du(1:nu)];
% general constraints
% lg = zeros(ng,1);
% ug = zeros(ng,1);
%if(ng>0)
%	if(nb==0)
%		for ii=1:min(nu,ng)
%			lg(ii) = -2.5; % lower bound
%			ug(ii) = -0.1; % - upper bound
%		end
%		for ii=min(nu,ng)+1:ng
%			lg(ii) = -10;
%			ug(ii) =  10;
%		end
%	else
%		for ii=1:ng
%			lg(ii) = -10;
%			ug(ii) =  10;
%		end
%	end
%end
%db(2*nu+1:end) = -4;
% llg = repmat(lg, 1, N);
% uug = repmat(ug, 1, N);
llg = reshape(lc,[ng,N]);
uug = reshape(uc,[ng,N]);
% general constraints last stage
lgN = zeros(ngN,1);
ugN = zeros(ngN,1);
if(ngN>0)
	if(nb==0)
%		for ii=1:min(nu,ng)
%			lg(ii) = -2.5; % lower bound
%			ug(ii) = -0.1; % - upper bound
%		end
		for ii=min(nu,ng)+1:ng
			lgN(ii) =   0;
			ugN(ii) =   0;
		end
	else
		for ii=1:ngN
			lgN(ii) =   0;
			ugN(ii) =   0;
		end
	end
end
%if(ng>0)
%	if(nb==0)
%		uub(1:min(nu,ng),N+1) =  0.0;
%	end
%end

% initial guess for states and inputs
x = zeros(nx, N+1); x(:,1) = ds0; % initial condition
u = zeros(nu, N);
mult_pi = zeros(nx,N);
mult_lam = zeros(2*(nb+ng)*N+2*(nb+ngN),1);
% mult_t = zeros(2*(nb+ng)*N+2*(nb+ngN),1);

free_x0 = 1; % consider x0 as optimization variable
warm_start = 0; % read initial guess from x and u

mu0 = 2;        % max element in cost function as estimate of max multiplier
kk = -1;		% actual number of performed iterations
k_max = 20;		% maximim number of iterations
tol = 1e-8;		% tolerance in the duality measure
infos = zeros(5, k_max);
inf_norm_res = zeros(1, 4);

tic;
HPIPM_d_solve_ipm2_hard_ocp_qp(kk, k_max, mu0, tol, N, nx, nu, nb, ng, ngN, time_invariant, free_x0, warm_start, AA, BB, bb, QQ, Qf, RR, SS, qq, qf, rr, llb, uub, CC, DD, llg, uug, CN, lgN, ugN, x, u, mult_pi, mult_lam, inf_norm_res, infos);
cpt_qp = toc*1e3;

dz = [x(:,1:N);u];
dxN = x(:,N+1);

lambda = mult_pi;

tmp = 2*(N+1)*nu;
mu = reshape(mult_lam(tmp+1:tmp+N*ng),[ng, N]) - reshape(mult_lam(tmp+N*ng+ngN+1:tmp+N*ng+ngN+N*ng),[ng, N]);

muN = mult_lam(tmp+N*ng+1:tmp+N*ng+ngN) - mult_lam(tmp+N*ng+ngN+N*ng:tmp+N*ng+ngN+N*ng+ngN);
end