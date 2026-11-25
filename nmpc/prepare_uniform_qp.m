function mem = prepare_uniform_qp(input, settings, mem)
%PREPARE_UNIFORM_QP Assemble the default QP using explicitly owned arrays.
if mem.sim_method == 1 && mem.hessian == 0 && mem.num_stages == 4
    qp = qp_generation(input, settings, mem);
    names = fieldnames(qp);
    for i = 1:numel(names)
        name = names{i};
        mem.(name) = qp.(name);
    end
else
    qp_generation(input, settings, mem);
end
end
