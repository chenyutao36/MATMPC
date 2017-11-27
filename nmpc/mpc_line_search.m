function [z_new, xN_new, lambda, mu, muN, mu_merit, alpha] = mpc_line_search(z,xN,dz,dxN,dw,od,y,yN,g,H,G,dG,B,dB,Q,QN,lambda_now,mu_now,muN_now, input,sizes,mpc_callnum)

    % for details refer to <<Numerical  Optimization>>, J. Nocedal and S. Wright, 2006 

    eta      = input.opt.meritfun.eta;
    tau      = input.opt.meritfun.tau;
    rho      = input.opt.meritfun.rho;
    mu_merit = input.opt.meritfun.mu_merit;
    mu_safty = input.opt.meritfun.mu_safty;
    
    alpha=1;
    if ~strcmp(input.opt.meritfun.second_order_correction,'yes') && ~strcmp(input.opt.meritfun.second_order_correction,'no')
        error('second_order_correction should be either yes or no');
    end

    if mpc_callnum>1
        l = norm(G,1)+sum( max(B,0) );   % constraint residual
        pd=dw'*H*dw;  % check if the Hessian along the specific direction is positive definite
        if pd>0
            sigma=0.5;
        else
            sigma=0;
        end
        mu_lb=(g'*dw+sigma*pd)/(1-rho)/l;
        if mu_merit<mu_lb
            mu_merit=mu_lb*mu_safty;
        end
        phi = GetObj( z,xN,od,y,yN,Q,QN,input,sizes )+mu_merit*l;
        D=g'*dw-mu_merit*l;
        newpoint=0;
        while ~newpoint && alpha>1e-4
            z_new   = z+alpha*dz;
            xN_new   = xN+alpha*dxN;
            [ G_new,B_new ] = GetGB( z_new, xN_new, od, input, sizes );
            l_new   =  norm( G_new,1 ) +sum( max(B_new, 0) );     
            phi_new = GetObj( z_new, xN_new, od, y, yN, Q, QN, input, sizes )+mu_merit*l_new;
            if phi_new<=phi+eta*alpha*D  % check if the merit function value is decreasing
                 newpoint=1;
            else 
                if strcmp(input.opt.meritfun.second_order_correction,'yes') && alpha==1                              
                    [ G_new,B_new ] = GetGB( z_new, xN_new, od, input, sizes );
                    deq   = G_new-dG*dw;
                    dineq = B_new-dB*dw;
                    input.opt.ipopt.options.ru=[-dineq;-deq];
                    input.opt.ipopt.options.rl=[-Inf*ones(sizes.nineq,1);-deq];
                    [dw_hat,fval_hat,exitflag_hat,info_ipopt_hat] = opti_ipopt(ipopt,[]);
                    dw_new=dw+dw_hat;
                    dz_new=reshape(dw_new(1:N*(nx+nu)), [nx+nu N])';
                    dxN_new=dw_new(N*(nx+nu)+1:end,1);
                    z_new   = z+dz_new;
                    xN_new   = xN+dxN_new;
                    [ G_new,B_new ] = GetGB( z_new, xN_new, od, input, sizes );
                    l_new   =  norm( G_new,1 ) +sum( max(B_new, 0) );   
                    phi_new = GetObj( z_new, xN_new, od, y, yN, Q, QN, input, sizes )+mu_merit*l_new;
                    if phi_new<=phi+eta*D
                        newpoint=1;
                    else
                        alpha=alpha*tau;
                    end % if
                else
                    alpha=alpha*tau;  % back tracking
                end % else if
            end % phi_new
        end % while
    end %if mpc_callnum

    z_new   = z+alpha*dz;
    xN_new  = xN+alpha*dxN;
    lambda  = alpha*lambda_now+(1-alpha)*input.lambda;
    mu      = alpha*mu_now+(1-alpha)*input.mu;
    muN     = alpha*muN_now+(1-alpha)*input.muN;
end