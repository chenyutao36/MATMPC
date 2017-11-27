function [Qh,S,R,A,B,gx,gu,Cx,Cu,c,a, tINT, tSENS] = mpc_qp_generation_cell(z,xN,y,yN,od,lambda,mu,muN,Q,QN,input,sizes)
%     tINT=0;
    tSENS=0;


%     N=sizes.N;
%     nx=sizes.nx;
%     nu=sizes.nu;
%     np=sizes.np;
%     ny=sizes.ny;
%     nyN=sizes.nyN;
%     nc=sizes.nc*2;
%     ncN=sizes.ncN*2;
%     
%     lb=input.lb;
%     ub=input.ub;
%     lbN=input.lbN;
%     ubN=input.ubN;
    
%     x_next=[z(1:nx,2:end),xN];
%     lambda_next=lambda(:,2:end);
%     
%     gx=input.data.gx;   
%     gu=input.data.gu;
%     a=input.data.a; 
%     Cx=input.data.Cx;
%     Cu=input.data.Cu;
%     c=input.data.c;
%     
%     Qh=input.data.Q;
%     S=input.data.S;
%     R=input.data.R;
%     A=input.data.A;
%     B=input.data.B;
         
    
%{        
    for i=1:N
           
        zi = z(:,i);
        refi = y(:,i);
        
        if ~isempty(od)
            parai = od(:,i);
        else
            parai = 0;
        end
        
%         switch integrator
%             case 'ERK4'
                tint=tic;        
                a(:,i)=F('F',zi,parai)-x_next(:,i);
                tINT=tINT+toc(tint)*1e3;
                
                 tsens=tic;
                 [Ai,Bi]=D('D',zi,parai);                                
                 tSENS=tSENS+toc(tsens)*1e3;
                    
                
%             case 'IRK4'
%                 
%                 [ xk_end, tsimu, sens_info] = IRK_INT_M(sizes, zi(1:nx), zi(nx+1:nx+nu), parai );
%                 tINT=tINT+tsimu;
%                 a(:,i)=xk_end-x_next(:,i);
%                 
%                 [ Sens, tsens] = IRK_SENS_M(sizes, sens_info );
%                 tSENS=tSENS+tsens;
%                 dG=Sens;
%         end   
        A{i}=full(Ai);
        B{i}=full(Bi);
        
%         switch hessian                    
%             case 'gauss_newton'
%                 dJi=full(dobji_vec_fun('jacobian_obji_vec_fun_0_0',zi,refi,parai,Q));
                [Jxi, Jui] = Ji_fun('Ji_fun',zi,parai,refi,Q);
                Jxi = full(Jxi);
                Jui = full(Jui);
%             case 'exact'
%                 mui = mu(:,i);
%                 switch integrator  
%                     case 'ERK4'                                                                                    
%                         Hi=Hi_EX('Hi_EX',zi,refi,parai,Q,lambda_next(:,i),mui);
%                             
%                     case 'IRK4'                           
%                         [ d2g, ~] = IRK_HESS_M(sizes, sens_info, lambda_next(:,i));
%                         d2f=dJ2dz_fun('dJ2dz_fun',zi,refi,parai,Q);
%                         d2b=dB2dz_fun('dB2dz_fun',zi,mui);
%                         Hi=d2f+sparse(d2g)+d2b;                           
%                 end
%         end
        
        
        Qh{i}=Jxi'*Jxi;
        S{i}=Jxi'*Jui;
        R{i}=Jui'*Jui;
        
        [gxi, gui]=gi_fun('gi_fun',zi,parai,refi,Q);
        gx(:,i)=full(gxi);
        gu(:,i)=full(gui);
                        
        ci=full(ineq_fun('ineq_fun',zi(1:nx), zi(nx+1:nx+nu), parai));
        c((i-1)*nc+1:i*nc,1)=[ci-ub(:,i);lb(:,i)-ci];
        
        [Cxi, Cui]=Ci_fun('Ci_fun',zi);
        Cx{i}=full(Cxi);
        Cu{i}=full(Cui);
    end
    
    i=N+1;
    zN = xN;
    refN = yN;
    if ~isempty(od)
         paraN = od(:,i);
     else
         paraN = 0;
    end
    
%     switch hessian
%             case 'gauss_newton'
                JN=JN_fun('JN_fun',zN,paraN,refN,QN);   
                JN = full(JN);
%             case 'exact'      
%                 switch integrator
%                     case 'ERK4'                     
%                         Hi=HN_EX('HN_EX',zN,refN,paraN,QN,muN);                                   
%                     case 'IRK4'                       
%                         d2fN=dJN2dz_fun('dJN2dz_fun',zN,refN,paraN,QN);
%                         d2bN=dBN2dz_fun('dBN2dz_fun',zN,muN);
%                         Hi=d2fN+d2bN;
%                         
%                 end          
%     end
   
    Qh{i}=JN'*JN;
              
    gx(:,N+1) = full(gN_fun('gN_fun',zN,paraN,refN,QN));
          
    cN=full(ineqN_fun('ineqN_fun',zN,paraN));
    c(N*nc+1:N*nc+ncN)=[cN-ubN;lbN-cN];
    
    Cx{i}=full(CN_fun('CN_fun',zN));
    %}

    %%
    tint = tic;
    [Qh,S,R,A,B,Cx,Cu,gx,gu,c,a] = qp_generation_mex(input, sizes);
%     [Qh,S,R,A,B,Cx,Cu,gx,gu,c,a] = qp_generation_mex(z, xN, y, yN, od, Q, QN, lb, ub, lbN, ubN, nx, nu, np, ny, nyN, nc, ncN, N);
    tINT = toc(tint)*1000;                      

end