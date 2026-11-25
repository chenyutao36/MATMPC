function mem = prepare_shooting_qp(input, settings, mem)
%PREPARE_SHOOTING_QP Assemble the QP for the configured prediction grid.
if mem.sim_method == 1 && mem.hessian == 0 && mem.num_stages == 4
    qp = qp_generation(input, settings, mem);
    names = fieldnames(qp);
    for i = 1:numel(names)
        name = names{i};
        mem.(name) = qp.(name);
    end
elseif isfield(mem, 'nonuniform_grid') && mem.nonuniform_grid
    qp_generation_ngrid(input, settings, mem);
else
    qp_generation(input, settings, mem);
end
end
