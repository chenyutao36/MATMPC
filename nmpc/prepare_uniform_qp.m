function mem = prepare_uniform_qp(input, settings, mem)
%PREPARE_UNIFORM_QP Assemble the QP using the shared shooting entry point.
mem = prepare_shooting_qp(input, settings, mem);
end
