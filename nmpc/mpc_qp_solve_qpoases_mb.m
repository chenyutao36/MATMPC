function [cpt_qp, mem] = mpc_qp_solve_qpoases_mb(sizes,mem,opt)

    nu=sizes.nu;
    N=sizes.N;  
    nbx=sizes.nbx;
    nc=sizes.nc;
    ncN=sizes.ncN;
    
    stageLower = reshape(mem.lb_du, nu, N);
    stageUpper = reshape(mem.ub_du, nu, N);
    blockLower = zeros(nu, mem.r);
    blockUpper = zeros(nu, mem.r);
    for block = 1:mem.r
        stages = mem.index_T(block)+1:mem.index_T(block+1);
        blockLower(:, block) = max(stageLower(:, stages), [], 2);
        blockUpper(:, block) = min(stageUpper(:, stages), [], 2);
    end
    
    initializing = mem.warm_start == 0;
    if initializing
        [mem.warm_start,sol,~,exitflag,~,multiplier,auxOutput] = qpOASES_sequence('i',mem.Hc_r,mem.gc_r,[mem.Ccx_r;mem.Ccg_r],...
            blockLower(:),blockUpper(:),[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
    elseif mem.hot_start == 1
        [sol,~,exitflag,~,multiplier,auxOutput] = qpOASES_sequence('h',mem.warm_start,mem.gc_r,...
            blockLower(:),blockUpper(:),[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
    else
        [sol,~,exitflag,~,multiplier,auxOutput] = qpOASES_sequence('m',mem.warm_start,mem.Hc_r,mem.gc_r,...
            [mem.Ccx_r;mem.Ccg_r],blockLower(:),blockUpper(:),[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
    end
    
    if exitflag ~= 0
        if initializing
            qpOASES_sequence('c', mem.warm_start);
        end
        error('MATMPC:QPSolverFailed', 'Move-blocking qpOASES failed (exitflag %d).', exitflag);
    end

    mem.du(:) = mem.T * sol;

    inputCount = mem.r * nu;
    generalOffset = inputCount + N * nbx;
    mem.mu_x_new(:) = -multiplier(inputCount+(1:N*nbx));
    mem.mu_new(:) = -multiplier(generalOffset+(1:N*nc+ncN));
    blockMu = reshape(-multiplier(1:inputCount), nu, mem.r);
    stageMu = zeros(nu, N);
    for block = 1:mem.r
        stages = mem.index_T(block)+1:mem.index_T(block+1);
        mu = blockMu(:, block);
        atLower = bsxfun(@eq, stageLower(:, stages), blockLower(:, block));
        atUpper = bsxfun(@eq, stageUpper(:, stages), blockUpper(:, block));
        active = bsxfun(@and, mu < 0, atLower) | bsxfun(@and, mu > 0, atUpper);
        % The lifted stage multipliers must sum to the reduced block multiplier.
        stageMu(:, stages) = bsxfun(@times, active, mu ./ max(1, sum(active, 2)));
    end
    mem.mu_u_new(:) = stageMu(:);
    cpt_qp   = auxOutput.cpuTime*1e3;
                
    Recover(mem, sizes);
end
