function [solution,cpt_qp,input] = mpc_qp_solve_simp(H,g,dB,B,input)

%     switch(input.opt.qpsolver)
%         case 'ipopt'
%                                   
%             H=sparse(H);
%             input.opt.ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
%             input.opt.ipopt.funcs.hessianstructure = @() tril(H);
%             input.opt.ipopt.funcs.objective=@(x) 0.5*x'*H*x + g'*x;
%             input.opt.ipopt.funcs.gradient = @(x) (H*x + g)';
%             input.opt.ipopt.options.A=sparse([dB;dG]);
%             input.opt.ipopt.options.ru=-[B;G];
%             input.opt.ipopt.options.rl=[-Inf*ones(nineq,1);-G];
%             [solution,fval,exitflag,info_ipopt] = opti_ipopt(input.opt.ipopt,[]);
%             lambda   = info_ipopt.Lambda.eqlin;
%             mu       = info_ipopt.Lambda.ineqlin;
%             info.cpt_qp   = info_ipopt.Time*1e3;
%             info.exitflag = exitflag;
                        
%         case 'qpoases'
            if input.opt.condensing_matrix.QP==0 
                [QP,solution,~,~,~,~,auxOutput] = qpOASES_sequence('i',H,g,dB,[],[],[],-B); 
                input.opt.condensing_matrix.QP=QP;
            else
                QP=input.opt.condensing_matrix.QP;
%                 if strcmp(input.opt.hotstart,'no')
                [solution,~,~,~,~,auxOutput] = qpOASES_sequence('m',QP,H,g,dB,[],[],[],-B);
%                 end
%                 if strcmp(input.opt.hotstart,'yes')
%                     [solution,~,~,~,~,auxOutput] = qpOASES_sequence('h',QP,g,[],[],[],-B);
%                 end
            end           
%             mu   = - multiplier(N*nu+1:end);
%             lambda = 0;
            cpt_qp   = auxOutput.cpuTime*1e3;
%             info.exitflag = exitflag;
            
%         case 'qpdunes'
%             H_dunes=cell2mat(H(1:N));
%             C=cell2mat(dG);
%             c=G;
%             
%             for i=1:numBlock
%                 zUpp = B(1:N*nc/2,1);
%                 zLow = B(N*nc/2+1:end,1);
%             end            
%             
%             if opt.iter==1
%                 qpDUNES( 'init', N, H_dunes, H{N+1}, g, C, c, zLow, zUpp, [], [], [], qpOptions );
%             else
%                 qpDUNES( 'update', H_dunes, H{N+1}, g, C, c, zLow, zUpp, [], [], [] );
%             end
%             [solution, stat, lambda, mu, fval] = qpDUNES( 'solve' );
%         case 'forces'
%             problem=sizes.problem;
%             
%             dB_full=full(dB(1:nineq,:));
%             H_full=full(H);
%             dG_full=full(dB(nineq+1:end,:));
%             G=B(nineq+1:end,1);
%             B=B(1:nineq,1);
%             
%             for i=1:N
%             
%                 a=G(i*nx+1:(i+1)*nx,1);
%                 Gxk=dG_full(i*nx+1:(i+1)*nx,(i-1)*nx+1:i*nx);
%                 Guk=dG_full(i*nx+1:(i+1)*nx,neq+(i-1)*nu+1:neq+i*nu);   
%                 gxk=g((i-1)*nx+1:i*nx,1);
%                 guk=g(neq+(i-1)*nu+1:neq+i*nu,1);
%                 Hk=blkdiag(H_full((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx), H_full(neq+(i-1)*nu+1:neq+i*nu,neq+(i-1)*nu+1:neq+i*nu) );
%                 
%                 problem=setfield(problem,['H',num2str(i)],Hk);
%                 problem=setfield(problem,['C',num2str(i)],[Gxk,Guk]);
%                 problem=setfield(problem,['a',num2str(i)],-a);
%                 problem=setfield(problem,['f',num2str(i)],[gxk;guk]);
%                 
%                 Ax=[dB_full((i-1)*nc+1:i*nc,i*nx+1:(i+1)*nx);dB_full(nc*N+(i-1)*nc+1:nc*N+i*nc,i*nx+1:(i+1)*nx)];
%                                
%                 b=-[B((i-1)*nc+1:i*nc,1);B(nc*N+(i-1)*nc+1:nc*N+i*nc,1)];
%                 
%                 if i<N
%                     Au=[dB_full((i-1)*nc+1:i*nc,neq+(i-1)*nu+1:neq+i*nu);dB_full(nc*N+(i-1)*nc+1:nc*N+i*nc,neq+(i-1)*nu+1:neq+i*nu)];
%                     problem=setfield(problem,['A',num2str(i+1)],[Ax,Au]);
%                     problem=setfield(problem,['b',num2str(i+1)],b); 
%                 else                    
%                     problem=setfield(problem,['A',num2str(i+1)],Ax);
%                     problem=setfield(problem,['b',num2str(i+1)],b); 
%                 end
%                 
%             end
% 
%             i=N+1;
%             Hk=H_full((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx);
%             gxk=g((i-1)*nx+1:i*nx,1);
%             problem=setfield(problem,['H',num2str(i)],Hk);
%             problem=setfield(problem,['f',num2str(i)],gxk);
%             
%             problem.ds0=G(1:nx,1);
%             [solverout,exitflag,info_forces] = MCA_ForcesPro(problem);
%             info.cpt_qp = info_forces.solvetime*1e3;
%             
%             triple=cell2mat(struct2cell(solverout));
%             solution=triple(1:nw,1);
%             info.lambda=triple(nw+1:nw+neq,1);
%             info.mu_now=triple(nw+neq+1:end,1);
%             info.exitflag = exitflag;
%             fval=info_forces.pobj;
%     end
end
