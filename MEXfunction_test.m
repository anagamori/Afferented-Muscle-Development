% Create a MEX configuration object
cfg = coder.config('mex');
% Turn on dynamic memory allocation
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% Generate MEX function
codegen -config cfg gridSpeedUpFcn -args {10}
disp('MEX generation complete!')

MToc = zeros(1,30);
MexToc = zeros(1,30);
for j = 1:30
    tic; gridSpeedUpFcn(10);
    MToc(j) = toc;
end

for j = 1:30
    tic; gridSpeedUpFcn_mex(10);
    MexToc(j) = toc;
end

eff_M_time = mean(MToc(1,21:30));
eff_Mex_time = mean(MexToc(1,21:30));

disp(['MATLAB code execution time (with nx=10): ', num2str(eff_M_time)])
disp(['MEX code executing time (with nx=10): ', num2str(eff_Mex_time)])

disp([' Speed up factor: ', num2str(eff_M_time/eff_Mex_time)]);