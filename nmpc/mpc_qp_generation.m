function [H,g,G,dB,B,dG,tINT,tSENS] = mpc_qp_generation(z,xN,y,yN,od,lambda,mu,muN,Q,QN,input,sizes)
    tINT=0;
    tSENS=0;


    N=sizes.N;
    nx=sizes.nx;
    nu=sizes.nu;
    
    x_next=[z(2:end,1:nx);xN'];
    lambda_next=lambda(2:end,:);
    
    H=input.data.H;    
    g=input.data.g;    
    G=input.data.G; 
    dG=input.data.dG;    
    B=input.data.B;
    dB=input.data.dB;
    
    lb=input.lb;
    ub=input.ub;
    lbN=input.lbN;
    ubN=input.ubN;
    
    integrator=input.opt.integrator;
    hessian=input.opt.hessian;  
    
    for i=1:N
           
        zi = z(i,:)';
        refi = y(i,:)';
        
        if ~isempty(od)
            parai = od(i,:)';
        else
            parai = 0;
        end
        
        switch integrator
            case 'ERK4'
                tint=tic;        
                G{i}=full(F('F',zi,parai))-x_next(i,:)';
                tINT=tINT+toc(tint)*1e3;
                
                tsens=tic;
                dG{i}=D('jacobian_F_0_0',zi,parai);
                tSENS=tSENS+toc(tsens)*1e3;
            case 'IRK4'
                
                [ xk_end, tsimu, sens_info] = IRK_INT_M(sizes, zi(1:nx), zi(nx+1:nx+nu), parai );
                tINT=tINT+tsimu;
                G{i}=xk_end-x_next(i,:)';
                
                [ Sens, tsens] = IRK_SENS_M(sizes, sens_info );
                tSENS=tSENS+tsens;
                dG{i}= sparse(Sens);
        end        
        
        switch hessian                    
                case 'gauss_newton'
                    dJi=dobji_vec_fun('jacobian_obji_vec_fun_0_0',zi,refi,parai,Q);
                    H{i}=dJi'*dJi;
                    
                case 'exact'
                    mui = mu(i,:)';
                    switch integrator  
                        case 'ERK4'                                                                                    
                            H{i}=Hi_EX('Hi_EX',zi,refi,parai,Q,lambda_next(i,:)',mui);
                            
                        case 'IRK4'                           
                            [ d2g, ~] = IRK_HESS_M(sizes, sens_info, lambda_next(i,:)' );
                            d2f=dJ2dz_fun('dJ2dz_fun',zi,refi,parai,Q);
                            d2b=dB2dz_fun('dB2dz_fun',zi,mui);
                            H{i}=d2f+sparse(d2g)+d2b;                           
                    end
        end
        
        g{i}=dJdz('dJdz',zi,refi,parai,Q);
                        
        bi=full(ineq_fun('ineq_fun',zi(1:nx), zi(nx+1:nx+nu), parai));
        B{i}=[bi;-bi]+[-ub(:,i);lb(:,i)];
        
        dB{i}=dBdz('dBdz',zi); 
    end
    
    i=N+1;
    zN = xN;
    refN = yN;
    if ~isempty(od)
         paraN = od(i,:)';
     else
         paraN = 0;
    end
    
    switch hessian
            case 'gauss_newton'
                dJN=dobjN_vec_fun('jacobian_objN_vec_fun_0_0',zN,refN,paraN,QN);
                H{i}=dJN'*dJN;                
            case 'exact'      
                switch integrator
                    case 'ERK4'                     
                        H{i}=HN_EX('HN_EX',zN,refN,paraN,QN,muN);                                   
                    case 'IRK4'                       
                        d2fN=dJN2dz_fun('dJN2dz_fun',zN,refN,paraN,QN);
                        d2bN=dBN2dz_fun('dBN2dz_fun',zN,muN);
                        H{i}=d2fN+d2bN;
                        
                end          
    end
              
    g{i}=dJNdz('dJNdz',zN,refN,paraN,QN);
          
    bN=full(ineqN_fun('ineqN_fun',zN,paraN));
    B{i}=[bN;-bN]+[-ubN;lbN];

    dB{i}=dBNdz('dBNdz',zN);
                  
end