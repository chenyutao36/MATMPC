classdef QpOasesMatrixCache < handle
    % Shared by MATLAB values that refer to the same native solver sequence.
    properties (SetAccess = private)
        Sequence
    end
    properties
        Hessian
        Constraints
        Valid = false
        Eligible = true
    end

    methods
        function obj = QpOasesMatrixCache(sequence)
            obj.Sequence = sequence;
        end
    end
end
