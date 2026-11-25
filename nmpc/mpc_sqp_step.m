function [input, mem, eq_res, ineq_res, KKT, OBJ] = mpc_sqp_step(input, settings, mem)
%MPC_SQP_STEP Select a step and update the iterate after acceptance.
[mem.alpha, mem.mu_merit] = Line_search(mem, input, settings);

if mem.alpha > 0
    input.x = input.x + mem.alpha * reshape(mem.dx, size(input.x));
    input.u = input.u + mem.alpha * reshape(mem.du, size(input.u));
    delta_lambda = mem.alpha * (mem.lambda_new(:) - input.lambda(:));
    mem.q_dual = reshape(delta_lambda, size(mem.q_dual));
    input.lambda = input.lambda + reshape(delta_lambda, size(input.lambda));

    [input.mu_u, delta_u] = updateMultipliers(input.mu_u, mem.mu_u_new, mem.alpha);
    [input.mu_x, delta_x] = updateMultipliers(input.mu_x, mem.mu_x_new, mem.alpha);
    [input.mu, delta_g] = updateMultipliers(input.mu, mem.mu_new, mem.alpha);
    mem.dmu = reshape([delta_u; delta_x; delta_g], size(mem.dmu));
end

if mem.alpha > 0 && settings.nz > 0
    [eq_res, ineq_res, KKT, OBJ, input.z] = solution_info(input, settings, mem);
else
    [eq_res, ineq_res, KKT, OBJ] = solution_info(input, settings, mem);
end
end

function [current, delta] = updateMultipliers(current, candidate, alpha)
% A Simulink bus can contain a scalar placeholder for an empty constraint set.
delta = candidate(:);
if isempty(delta)
    return;
end
delta = alpha * (delta - current(:));
current = current + reshape(delta, size(current));
end
