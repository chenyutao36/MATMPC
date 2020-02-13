classdef qpalm < handle
    % qpalm interface class for QPALM solver 
    % This class provides a complete interface to the C implementation
    % of the QPALM solver.
    %
    % qpalm Methods:
    %
    %   setup             - configure solver with problem data
    %   solve             - solve the QP
    %   warm_start        - set warm starting variables x and y 
    %
    %   default_settings  - create default settings structure
    %   current_settings  - get the current solver settings structure

    properties (SetAccess = private, Hidden = true)
        n;
        m;
%         objectHandle % Handle to underlying C instance
    end
   methods(Static) 
        %%
        function out = default_settings()
            % DEFAULT_SETTINGS get the default solver settings structure
            out = qpalm_mex('default_settings');
        end
        
        %%
        function out = constant(constant_name)
            % CONSTANT Return solver constant
            %   C = CONSTANT(CONSTANT_NAME) return constant called CONSTANT_NAME
            out = qpalm_mex('constant', constant_name);
        end
        
    end
    methods
        %% Constructor - Create a new solver instance
        function this = qpalm(varargin)
            % Construct QPALM solver class
%             this.objectHandle = qpalm_mex('new', varargin{:});
        end

        %% Destructor - destroy the solver instance
        function delete(this)
            % Destroy QPALM solver class
            qpalm_mex('delete');
        end
            
            
        %% Setup function
        function varargout = setup(this, varargin)
            % SETUP configure solver with problem data
            %
            %   setup(Q,q,A,bmin,bmax,options) or
            %   setup(Q,q,A,bmin,bmax,x,y,options)

            nargin = length(varargin);

            %dimension checks on user data. Mex function does not
            %perform any checks on inputs, so check everything here
            assert(nargin >= 5, 'incorrect number of inputs');
            assert(nargin ~= 7, 'incorrect number of inputs');
            [Q,q,A,bmin,bmax] = deal(varargin{1:5});
            
            %
            % Get problem dimensions
            %

            % Get number of variables n
            if (isempty(Q))
                if (~isempty(q))
                    n = length(q);
                else
                    if (~isempty(A))
                        n = size(A, 2);
                    else
                        error('The problem does not have any variables');
                    end
                end
            else
                n = size(Q, 1);
                assert(n==size(Q,2), 'Q must be a square matrix');
            end
            
            this.n = n;

            % Get number of constraints m
            if (isempty(A))
                m = 0;
            else
                m = size(A, 1);
            end

            this.m = m;
            %
            % Create sparse matrices and full vectors if they are empty
            %

            if (isempty(Q))
                Q = sparse(n, n);
%                 Q = zeros(n,n);
            else
                Q   = sparse(Q);
%                 Q = full(Q(:,:));
            end
            if (isempty(q))
                q = zeros(n, 1);
            else
                q   = full(q(:));
            end

            % Create proper constraints if they are not passed
            if (isempty(A) && (~isempty(bmin) || ~isempty(bmax))) || ...
                (~isempty(A) && (isempty(bmin) && isempty(bmax)))
                error('A must be supplied together with at least one bound l or u');
            end

            if (~isempty(A) && isempty(bmin))
                bmin = -Inf(m, 1);
            end

            if (~isempty(A) && isempty(bmax))
                bmax = Inf(m, 1);
            end

            if (isempty(A))
                A = sparse(m, n);
%                 A = zeros(m,n);
                bmin = -Inf(m, 1);
                bmax = Inf(m, 1);
            else
                bmin = full(bmin(:));
                bmax = full(bmax(:));
                A = sparse(A);
%                 A = full(A(:,:));
            end


            %
            % Check vector dimensions (not checked from the C solver)
            %

            assert(length(q)    == n, 'Incorrect dimension of q');
            assert(length(bmin) == m, 'Incorrect dimension of l');
            assert(length(bmax) == m, 'Incorrect dimension of u');

            %
            % Convert infinity values to QPALM_INFINITY
            %
            bmax = min(bmax, qpalm.constant('QPALM_INFTY'));
            bmin = max(bmin, -qpalm.constant('QPALM_INFTY'));


            %make a settings structure from the remainder of the arguments.
            %'true' means that this is a settings initialization, so all
            %parameter/values are allowed.  No extra inputs will result
            %in default settings being passed back
            if nargin == 5
                theSettings = validateSettings([]);
            elseif nargin == 6
                theSettings = validateSettings(varargin{6:end});
            elseif nargin >= 8
                theSettings = validateSettings(varargin{8:end});
                theSettings.warm_start = true;
            end

            [varargout{1:nargout}] = qpalm_mex('setup',n,m,Q,q,A,bmin,bmax,theSettings);
            
            if nargin >= 8
                [x,y] = deal(varargin{6:7});
                this.warm_start(x,y);
            end

        end

        function [n,m] = get_dimensions(this)
           n = this.n;
           m = this.m;
           return;
        end

        %% Warm start primal and dual variables
        function warm_start(this, x, y)
            % WARM_START warm start primal and/or dual variables
            %
            %   warm_start(x,y)

            % Get problem dimensions
            [n, m]  = get_dimensions(this);
            assert(isempty(x) || length(x) == n, 'Incorrect dimension of x');
            assert(isempty(y) || length(y) == m, 'Incorrect dimension of y');

            if (~isempty(x))
                x = full(x(:));
            end
            if (~isempty(y))
                y = full(y(:));
            end

            qpalm_mex('warm_start',x,y);

        end
        
        %% Update functions
        function update_q(this, q)
            % UPDATE_Q update the linear part of the cost
            %
            %   update_q(q)
            
            if (~isempty(q))
                q = full(q(:));
            end
            qpalm_mex('update_q', q);
        end
        
        function update_bounds(this, bmin, bmax)
            % UPDATE_BOUNDS update the lower and upper bounds of the linear
            % constraints. Use empty matrices to denote bounds that are not
            % updated.
            %
            %   update_bounds(bmin, bmax)
            
            if (~isempty(bmin))
                bmin = full(bmin(:));
            end
            if (~isempty(bmax))
                bmax = full(bmax(:));
            end
            bmax = min(bmax, qpalm.constant('QPALM_INFTY'));
            bmin = max(bmin, -qpalm.constant('QPALM_INFTY'));
            
            qpalm_mex('update_bounds', bmin, bmax);
        end
            

        %% Solve function
        function varargout = solve(this, varargin)
            % SOLVE solve the QP

            nargoutchk(0,1);  %either return nothing (but still solve), or a single output structure
            [out.x, out.y, out.prim_inf_cert, out.dual_inf_cert, out.info] = qpalm_mex('solve');
            if(nargout)
                varargout{1} = out;
            end
%             delete(this);
            return;
        end

    end
end



function currentSettings = validateSettings(varargin)

%get the default settings
currentSettings = qpalm_mex('default_settings');

%no settings passed -> return defaults
if(isempty(varargin))
    return;
end

%check for structure style input
if(isstruct(varargin{1}))
    newSettings = varargin{1};
    assert(length(varargin) == 1, 'too many input arguments');
else
    newSettings = struct(varargin{:});
end

%get the qpalm settings fields
currentFields = fieldnames(currentSettings);

%get the requested fields in the update
newFields = fieldnames(newSettings);

%check for unknown parameters
badFieldsIdx = find(~ismember(newFields,currentFields));
if(~isempty(badFieldsIdx))
    error('Unrecognized solver setting ''%s'' detected',newFields{badFieldsIdx(1)});
end



%everything checks out - merge the newSettings into the current ones
for i = 1:length(newFields)
    currentSettings.(newFields{i}) = double(newSettings.(newFields{i}));
end


end


