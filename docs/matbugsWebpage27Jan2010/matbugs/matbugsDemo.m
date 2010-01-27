% An example illustrating how to call matbugs.m
% We use the schools data/model from p140/ p592 of Gelman's book
% Written by Maryam Mahdaviani, August 2005
% Modified by Kevin Murphy (murphyk@cs.ubc.ca), 10 October 2005

dataStruct = struct('J', 8, ...
                    'y', [28, 8, -3, 7, -1, 1, 18, 12], ...
                   'sigma_y', [15, 10, 16, 11, 9, 11, 10, 18]);

Nchains = 3;

% we initialize the params to the observed data values, but with decreasing confidence,
% as suggested on p593 of Gelman
for i=1:Nchains
  S.theta = dataStruct.y;
  S.mu_theta = 0;
  S.sigma_theta = 10^i; % each chain becomes more over-dispersed
  initStructs(i) = S;
end



[samples, stats] = matbugs(dataStruct, ...
		fullfile(pwd, 'schools_model.txt'), ...
		'init', initStructs, ...
		'view', 0, 'nburnin', 1000, 'nsamples', 500, ...
		'thin', 10, ...
		'monitorParams', {'theta', 'mu_theta', 'sigma_theta'}, ...
		'Bugdir', 'C:/Program Files/WinBUGS14');

disp('Rhat')
stats.Rhat

disp('means')
stats.mean

figure;
clf
N = 8+2; % monitor 10 variables
colors = 'rgb';
for j=1:8
  subplot(3,3,j); hold on
  for c=1:Nchains
    plot(samples.theta(c,:,j), colors(c));
  end
  title(sprintf('theta %d', j));
end
subplot(3,3,9); hold on
for c=1:Nchains
  plot(samples.mu_theta(c,:), colors(c));
end
title(sprintf('mu.theta'))
