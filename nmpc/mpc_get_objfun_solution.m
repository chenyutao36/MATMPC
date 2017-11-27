
function [KKT,eq_res,ineq_res] = mpc_get_objfun_solution(lambda,mu,muN, ds0, input,sizes) 
    
    %{
    Q = input.W;
    QN = input.WN;

    N=sizes.N;
    nx=sizes.nx;
    nu=sizes.nu;
    np=sizes.np;
    ny=sizes.ny;
    nyN=sizes.nyN;
    nc=2*sizes.nc;
    ncN=2*sizes.ncN;
    
    lb=input.lb;
    ub=input.ub;
    lbN=input.lbN;
    ubN=input.ubN;
       
    gx=input.data.gx;   
    gu=input.data.gu;
    a=input.data.a; 
    Cx=input.data.Cx;
    Cu=input.data.Cu;
    c=input.data.c;
    dL = input.data.dL;   
    x_next=[z(1:nx,2:end),xN];
          
    for i=1:N
           
        zi = z(:,i);
        refi = y(:,i);
        
        if ~isempty(od)
            parai = od(:,i);
        else
            parai = 0;
        end
        
        lambdai = lambda(:,i);              
        mui = mu(:,i);
        lambda_next = lambda(:,i+1);
              
        a(:,i)=F('F',zi,parai)-x_next(:,i);
        
        ci=full(ineq_fun('ineq_fun',zi(1:nx), zi(nx+1:nx+nu), parai));
        c((i-1)*nc+1:i*nc,1)=[ci-ub(:,i);lb(:,i)-ci];    
        
        [dobji, adj_dGi, adj_dBi] = adj_fun('adj_fun',zi,parai,refi,Q, lambda_next, mui);
        dL{i} = dobji + adj_dGi-[lambdai;zeros(nu,1)]+adj_dBi;

    end
    
    i=N+1;
    zN = xN;
    refN = yN;
    if ~isempty(od)
         paraN = od(:,i);
     else
         paraN = 0;
    end
         
    cN=full(ineqN_fun('ineqN_fun',zN,paraN));
    c(N*nc+1:N*nc+ncN)=[cN-ubN;lbN-cN];
    
    [dobjN, adj_dBN] = adjN_fun('adjN_fun',zN,paraN,refN,QN, muN);
    dL{i} = dobjN +adj_dBN;
    
    %%
    dL=full(cell2mat(dL));
    
    eq_res=norm( [ds0;reshape(a,[N*nx,1])],1 );
    ineq_res=sum( max( c,0 ) );
    KKT=norm( dL,1 ); 
    
%     Termination=max([KKT, eq_res, ineq_res]);
        
%}
%%

%     [eq_res, ineq_res, KKT] = solution_info(z, xN, y, yN, od, Q, QN, lb, ub, lbN, ubN, lambda, mu, muN, ds0, nx, nu, np, ny, nyN, nc, ncN, N);
    [eq_res, ineq_res, KKT] = solution_info(lambda, mu, muN, ds0, input, sizes);
end