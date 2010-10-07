bugsFolder = 'C:\kmurphy\Programs\WinBUGS14';
%demoFolder = 'C:\kmurphy\Programs\WinBUGS14\Manuals\Tutorial';
demoFolder = 'C:\kmurphy\GoogleCode\matbugs\demos\seedsSimpleDemo';

% read data
% r n x1 x2
data = importdata(fullfile(demoFolder, 'seeds_data.txt'), '\t', 1);
data.data(:,1) = []; % strip off empty column
dataStruct = struct('N', size(data.data,1), ...
  'r', data.data(:,1), ...
  'n', data.data(:,2), ...
  'x1', data.data(:,3), ...
  'x2', data.data(:,4));
 
Nchains = 2;
% specify initial values by cutting and pasting from init file
S.alpha0 = 0; S.alpha1 = 0; S.alpha2 = 0; S.alpha12 = 0;
S.tau = 1;
%S.sigma = 2;
S.b = zeros(1,21); 
initStructs(1) = S;

S.alpha0 = 10; S.alpha1 = 10; S.alpha2 = 10; S.alpha12 = 10;
S.tau = 0.5;  
%S.sigma= 5;
S.b = [0.1, -0.2, 0.25, 0.11, -0.21, 0.3, -0.25, 0.15, -0.31, -0.1,...
0.1, 0.12, 0.2, -0.2, 0.4, -0.24, 0.14, 0.3, -0.2, 0.1, 0.05];
initStructs(2) = S;

% run sampler
[samples, stats] = matbugs(dataStruct, ...
  fullfile(demoFolder, 'seeds_model.txt'), ...
  'init', initStructs, ...
  'view', 1, 'nburnin', 1000, 'nsamples', 500, ...
  'thin', 10, ...
  'monitorParams', {'alpha0','alpha1','alpha12'}, ...
  'Bugdir', bugsFolder);
