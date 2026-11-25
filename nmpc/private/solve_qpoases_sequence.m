function [sol, multiplier, elapsed, mem] = solve_qpoases_sequence(mem, H, gradient, A, ...
    lower, upper, lowerA, upperA, solverName)
%SOLVE_QPOASES_SEQUENCE Choose a sequence call based on the installed QP matrices.
initializing = mem.warm_start == 0;
checkHotstart = mem.hot_start == 1;
cache = [];
shortHotstart = false;
if ~initializing && isfield(mem, 'qpoases_matrix_cache') && ...
        mem.qpoases_matrix_cache.Sequence == mem.warm_start
    cache = mem.qpoases_matrix_cache;
    shortHotstart = checkHotstart && cache.Eligible && cache.Valid && ...
        isequal(H, cache.Hessian) && isequal(A, cache.Constraints);
    % A failed native call may have changed the sequence even though MATLAB
    % still has the previously returned mem. Invalidate that shared record.
    cache.Valid = false;
end
if isfield(mem, 'qpoases_matrix_cache') && (~checkHotstart || (~isempty(cache) && ~cache.Eligible))
    if ~isempty(cache)
        cache.Eligible = false;
        cache.Hessian = [];
        cache.Constraints = [];
    end
    mem = rmfield(mem, 'qpoases_matrix_cache');
end

if initializing
    [mem.warm_start, sol, ~, exitflag, ~, multiplier, aux] = qpOASES_sequence('i', ...
        H, gradient, A, lower, upper, lowerA, upperA, mem.qpoases_opt);
elseif shortHotstart
    [sol, ~, exitflag, ~, multiplier, aux] = qpOASES_sequence('h', mem.warm_start, ...
        gradient, lower, upper, lowerA, upperA, mem.qpoases_opt);
else
    [sol, ~, exitflag, ~, multiplier, aux] = qpOASES_sequence('m', mem.warm_start, ...
        H, gradient, A, lower, upper, lowerA, upperA, mem.qpoases_opt);
end

if exitflag ~= 0
    if initializing
        qpOASES_sequence('c', mem.warm_start);
    end
    error('MATMPC:QPSolverFailed', '%s failed (exitflag %d).', solverName, exitflag);
end
if initializing && isfield(mem, 'own_qpoases_sequence') && mem.own_qpoases_sequence
    sequence = mem.warm_start;
    mem.qpoases_cleanup = onCleanup(@() qpOASES_sequence('c', sequence));
end

% Only sequences initialized with opt-in can acquire a cache. Otherwise a
% previous MATLAB value with no shared record could update the same handle.
if checkHotstart && (initializing || (~isempty(cache) && cache.Eligible))
    if isempty(cache)
        cache = matmpc_internal.QpOasesMatrixCache(mem.warm_start);
    end
    if ~shortHotstart
        cache.Hessian = copyMatrix(H);
        cache.Constraints = copyMatrix(A);
    end
    cache.Valid = true;
    mem.qpoases_matrix_cache = cache;
end
elapsed = aux.cpuTime * 1e3;
end

function copy = copyMatrix(matrix)
% A plain assignment can share data with the original MATLAB matrix.
% Padding before taking the slice also detaches scalar and vector storage.
if isempty(matrix)
    copy = matrix;
else
    padded = [matrix(:); matrix(1)];
    copy = reshape(padded(1:end-1), size(matrix));
end
end
