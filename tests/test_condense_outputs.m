function tests = test_condense_outputs
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler with mex -setup C to run the condensing gateway test.');
root = fileparts(fileparts(mfilename('fullpath')));
gateway = fullfile(root, 'mex_core');
core = fullfile(root, 'core');
build = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(build.Folder));
testCase.addTeardown(@unloadFixture);
mex('-silent', '-outdir', build.Folder, ['-I' core], '-output', 'explicit_condensing', ...
    fullfile(gateway, 'Condensing.c'), fullfile(gateway, 'matmpc_la_mex.c'), ...
    fullfile(core, 'condense_matrix.c'), fullfile(core, 'condense_affine.c'), '-lmwblas');
end

function unloadFixture
clear explicit_condensing
end

function testExplicitOutputWithoutWorkspaceFields(testCase)
settings = struct('nx',1,'nu',1,'N',2,'nbx',1,'nc',0,'ncN',1);
mem = struct('Q',ones(1,3),'S',zeros(1,2),'R',ones(1,2),'A',ones(1,2),'B',ones(1,2), ...
    'Cx',1,'Cgx',zeros(0,2),'Cgu',zeros(0,2),'CgN',2,'gx',zeros(1,3),'gu',[0.3,-0.2], ...
    'a',[0.25,-0.10],'ds0',0.40,'lc',-4,'uc',7,'lb_dx',[-3;-2],'ub_dx',[2;3]);
before = mem;
qp = explicit_condensing(mem, settings);
state_map = tril(ones(2));
offset = [0.65;0.55];
verifyEqual(testCase, mem, before);
verifyEqual(testCase, qp.Hc, eye(2) + state_map' * state_map, 'AbsTol', 1e-12);
verifyEqual(testCase, qp.gc, mem.gu' + state_map' * offset, 'AbsTol', 1e-12);
verifyEqual(testCase, qp.Ccx, state_map);
verifyEqual(testCase, qp.Ccg, [2,2]);
verifyEqual(testCase, [qp.lcc,qp.ucc], [mem.lc,mem.uc] - 2*offset(end), 'AbsTol', 1e-12);
verifyEqual(testCase, [qp.lxc,qp.uxc], [mem.lb_dx,mem.ub_dx] - repmat(offset,1,2), 'AbsTol', 1e-12);
end

function testRejectsInvalidDimensionsAndBuffers(testCase)
verifyError(testCase, @() explicit_condensing(), 'MATMPC:Condensing:Arguments');
settings = struct('nx',1,'nu',1,'N',2,'nbx',1,'nc',0,'ncN',0);
verifyError(testCase, @() explicit_condensing(struct(),settings), 'MATMPC:Condensing:Buffer');
settings.N = 1.5;
verifyError(testCase, @() explicit_condensing(struct(),settings), 'MATMPC:Condensing:Dimensions');
end
