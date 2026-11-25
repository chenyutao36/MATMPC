function smoke_simulink_inverted_pendulum()
%SMOKE_SIMULINK_INVERTED_PENDULUM Build and exercise both bundled charts.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
originalPath = path;
restorePath = onCleanup(@() path(originalPath));
addpath(projectRoot);
Build_Model('InvertedPendulum');
variables = initialVariables();

testController(projectRoot, 'MATMPC_SIMULINK', variables);
testController(projectRoot, 'MATMPC_SIMULINK_R2017b', variables);
end

function variables = initialVariables()
settings = struct();  % The script populates this instead of the MATLAB settings function.
Initialization_Simulink;
assert(strcmp(opt.RTI, 'yes'), 'The default Simulink controller should use RTI.');
variables = struct('settings', settings, 'opt', opt, 'Ts', Ts, 'x0', x0, 'nx', nx, ...
    'nu', nu, 'ny', ny, 'nyN', nyN, 'N', N, 'z0', z0, 'lb', lb, 'ub', ub, ...
    'lbx', lbx, 'ubx', ubx, 'lbu', lbu, 'ubu', ubu, 'W', W, 'WN', WN, 'para', para);
end

function testController(projectRoot, model, variables)
assert(~bdIsLoaded(model), 'Close %s before running this smoke test.', model);
load_system(fullfile(projectRoot, [model, '.slx']));
closeModel = onCleanup(@() close_system(model, 0));

controller = get_param([model '/MATLAB Function'], 'PortHandles');
signals = {'control', 'solveTime', 'residual', 'objective'};
for i = 1:numel(signals)
    set_param(controller.Outport(i), 'DataLogging', 'on', ...
        'DataLoggingNameMode', 'Custom', 'DataLoggingName', signals{i});
end

testInput = Simulink.SimulationInput(model);
names = fieldnames(variables);
for i = 1:numel(names)
    name = names{i};
    testInput = testInput.setVariable(name, variables.(name));
end
Ts = variables.settings.Ts_st;
testInput = testInput.setModelParameter('StopTime', num2str(2 * Ts), ...
    'FixedStep', num2str(Ts), 'SignalLogging', 'on');
result = sim(testInput);

for i = 1:numel(signals)
    samples = result.logsout.get(signals{i}).Values.Data;
    assert(numel(samples) >= 2 && all(isfinite(samples(:))), ...
        '%s must produce successive finite %s samples.', model, signals{i});
end
controls = result.logsout.get('control').Values.Data;
assert(all(controls(:) >= variables.lbu(1) - 1e-6 & ...
    controls(:) <= variables.ubu(1) + 1e-6), '%s control bound violated.', model);
fprintf('%s passed the Simulink controller smoke test.\n', model);
end
