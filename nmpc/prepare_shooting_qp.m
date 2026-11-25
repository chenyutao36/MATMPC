function mem = prepare_shooting_qp(input, settings, mem, moveBlocking)
%PREPARE_SHOOTING_QP Assemble the QP for the configured prediction grid and control blocks.
if nargin < 4
    moveBlocking = false;
end
if mem.sim_method == 1 && mem.hessian == 0 && mem.num_stages == 4
    if moveBlocking
        qp = qp_generation(input, settings, mem, mem.index_T);
    else
        qp = qp_generation(input, settings, mem);
    end
    names = fieldnames(qp);
    for i = 1:numel(names)
        name = names{i};
        mem.(name) = qp.(name);
    end
elseif moveBlocking
    qp_generation_mb(input, settings, mem);
elseif isfield(mem, 'nonuniform_grid') && mem.nonuniform_grid
    qp_generation_ngrid(input, settings, mem);
else
    qp_generation(input, settings, mem);
end
end
