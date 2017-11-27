function [ G,B ] = GetGB( z,xN,od,input,sizes )

    N=sizes.N;
    nx=sizes.nx;
    nu=sizes.nu;
        
    G=input.data.G; 
    B=input.data.B;
        
    x_next=[z(2:end,1:nx);xN'];
    
    lb=input.lb;
    ub=input.ub;
    lbN=input.lbN;
    ubN=input.ubN;
    
    integrator=input.opt.integrator; 
    
    for i=1:N
           
        zi = z(i,:)';
        if ~isempty(od)
            parai = od(i,:)';
        else
            parai = 0;
        end
        
        switch integrator
            case 'ERK4'
                G{i}=full(F('F',zi,parai))-x_next(i,:)';
            case 'IRK4'
                [ xk_end, ~, ~] = IRK_INT_M(sizes, zi(1:nx), zi(nx+1:nx+nu), parai );
                G{i}=xk_end-x_next(i,:)';
        end               
        
        bi=full(ineq_fun('ineq_fun',zi(1:nx), zi(nx+1:nx+nu), parai));
        B{i}=[bi;-bi]+[-ub(:,i);lb(:,i)];
        
    end
    
    i=N+1;
    zN = xN;
    
    if ~isempty(od)
         paraN = od(i,:)';
     else
         paraN = 0;
    end
    
    bN=full(ineqN_fun('ineqN_fun',zN,paraN));
    B{i}=[bN;-bN]+[-ubN;lbN];
      
    G=[input.x0-z(1,1:nx)';cell2mat(G)];
    B=cell2mat(B);

end

