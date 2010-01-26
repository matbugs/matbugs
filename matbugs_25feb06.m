function [samples, stats, structArray] = matbugs(dataStruct, bugsModel, varargin)
% MATBUGS a Matlab interface for WinBugs, similar to R2WinBUGS
% [samples, stats] = matbugs(dataStruct,  bugsModelFileName, ...)
%
% This generates a file called 'script.txt' in the directory where
% winbugs14.exe is kept (so you must have write access to this directory!),
% then calls winbugs with this script, then reads the resulting samples from
% files into matlab structs.
%
% INPUT:
% dataStruct contains values of observed variables.
% bugsModel is the name of the model file; MUST END IN .txt
% 
% Note: variables with names 'a.b' in the bugs model file
% should be called 'a_b' in the matlab data structure.
%
% Optional arguments passed as 'string', value pairs [default in brackets, case
% insensitive]:
%
% 'init' - init(i).v is a struct containing initial values for variable 'v'
%          for chain i.  Uninitialized variables are given random initial
%          values by WinBUGS.  It is highly recommended that you specify the
%          initial values of all root (founder) nodes.
% 'monitorParams' - cell array of field names (use 'a_b' instead of 'a.b')
%                   [defaults to *, which currently does nothing...]
% 'nChains'  - number of chains [3]
% 'nBurnin'  - num samples for burn-in per chain [1000]
% 'nSamples' - num samples to keep after burn-in [5000]
% 'thin'     - keep every n'th step [1]
% 'blocking' - do block updating for GLMs (using IRLS proposal)?
%              [1 - doesn't yet work]
% 'view'     - set to 1 if you want to view WinBUGS output (then close the WinBUGS
%                window to return control to matlab)
%              set to 0 to close WinBUGS automatically without pausing [default 0]
% 'Bugdir'   - location of winbugs executable ['C:/Program Files/WinBUGS14']
% 'workingDir' - directory to store temporary data/init/coda files [pwd/tmp]
%
% Note that total number of iterations  = nBurnin + nSamples * thin.
%
% OUTPUT
% S contains the samples; each field may have a different shape:
%  S.theta(c, s)       is the value of theta in sample s, chain c
%                      (scalar variable)
%  S.theta(c, s, i)    is the value of theta(i) in sample s, chain c
%                      (vector variable)
%  S.theta(c, s, i, j) is the value of theta(i,j) in sample s, chain c
%                      (matrix variable)
%
% stats contains various statistics, currently:
%    stats.mean, stats.std and stats.Rhat.
% Each field may have a different shape:
%    stats.mean.theta
%    stats.mean.theta(i)
%    state.mean.theta(i,j)
%
% Rhat is the "estimated potential scale reduction" statistic due to
%     Gelman and Rubin.
%
% Rhat values less than 1.1 mean the chain has probably converged for
%     this variable.
%
% Example
%
%     [S,stats] = matbugs(data, 'school_model.txt', 'nSamples', 10000, ..
%                          'monitorParams', {'mu_theta', 'theta'}, ...
%                          'init', initStruct, 'nchains', 3)
%
% Written by Maryam Mahdaviani (maryam@cs.ubc.ca)
% and Kevin Murphy (murphyk@cs.ubc.ca), August 2005
%
% Last updated 19 December 2005.


[initStructs, Bugdir, nChains, view, workingDir, nBurnin, nSamples, monitorParams, thin, blocking, junk] =  ...
   process_options(...
	varargin, ...
	'init', {}, ...
	'Bugdir', 'C:/Program Files/WinBUGS14', ...
	'nChains', 3, ...
	'view', 0, ...
	'workingDir', fullfile(pwd,'tmp'), ...
	'nBurnin', 1000, ...
	'nSamples', 5000, ...
	'monitorParams', {}, ...
	'thin', 1, ...
	'blocking', 1);

if length(initStructs) ~= nChains
  error(['init structure does not match number of chains ', ...
    sprintf('(%d)', nChains)]);
end

if ~exist(workingDir, 'dir')
  mkdir(pwd, 'tmp');
end

log_filename = fullfileKPM(workingDir, 'log.txt');
his_filename = fullfileKPM(workingDir, 'history.txt');

scriptFile = [Bugdir,'\','script.txt'];

%bugsModel = [pwd, '/', bugsModel];
%bugsModel = fullfile(pwd, bugsModel);
bugsModel = strrep(bugsModel, '\', '/'); % winBUGS wants a/b not a\b

codaFile = fullfileKPM(workingDir, 'coda');

fid = fopen(scriptFile,'w');
if (fid == -1)
   error(['Cannot open ', scriptFile]);
end
%display(option)
fprintf(fid, 'display(''log'') \n');

%check(model file)
fprintf(fid, 'check(''%s'')\n',bugsModel);

if ~isempty(dataStruct)
  dataFile = fullfileKPM(workingDir, 'data.txt');
  dataGen(dataStruct, dataFile);
  fprintf(fid, 'data(''%s'')\n', dataFile);
end
%fprintf(fid, '.txt'')\n');

fprintf(fid, 'compile(%u) \n', nChains);

initfileN = size(initStructs,2);
for i=1:initfileN
  initFileName = fullfileKPM(workingDir, ['init_', num2str(i) '.txt']);
  dataGen(initStructs(i), initFileName)
  fprintf(fid, 'inits (%u, ''%s'')\n', i, initFileName);
end

if 0 % blocking
  fprintf(fid, 'blockfe(1)\n'); % block fixed effects (for GLMs)
end
fprintf(fid, 'gen.inits() \n');
fprintf(fid, 'update(%u)\n', nBurnin);

%setting params to monitor
if isempty(monitorParams)
  fprintf(fid, 'set (*)\n');
else
  for i=1:length(monitorParams)
    fprintf(fid, 'set (%s)\n', strrep(monitorParams{i}, '_', '.'));
  end
  %fprintf(fid, 'set (%s)\n','deviance');
end


%fprintf(fid, 'dic.set()\n');
fprintf(fid, 'thin.updater(%u)\n', thin);
fprintf(fid, 'update(%u)\n', nSamples);
fprintf(fid, 'coda(*, ''%s'')\n',  codaFile);
fprintf(fid, 'stats(*)\n');


%fprintf(fid, 'dic.stats()\n');
fprintf(fid, 'history(*, ''%s'')\n', his_filename);
fprintf (fid, 'save (''%s'')\n',  log_filename);

if (view == 0)
  fprintf(fid, 'quit() \n');
end
fclose(fid);

%calling WinBUGS and passing the script to it
str = ['"', Bugdir, '\Winbugs14.exe" /PAR script.txt'];
dos(str); % DOS wants a\b

% passing the coda files to bugs2mat to convert files to structs
codaIndex = [codaFile, 'Index.txt'];
for i=1:nChains
  codaF = [codaFile, num2str(i), '.txt'];
  S = bugs2mat(codaIndex, codaF);
  structArray(i) = S;
end

samples = structsToArrays(structArray);
stats = computeStats(samples);

if nChains == 1
   disp('EPSR not calculated (only one chain)');
end

%%%%%%%%%%%%%

function dataGen(dataStruct, fileName)
% This is a helper function to generate data or init files for winBUGS
% Inputs:
%   fileName: name of the text file containing initial values. for each
%             chain we'll fileName_i where 'i' is the chain number,
%   dataStruct: is a Struct with name of params(consistant in the same
%               order with paramList) are fields and intial values are functions

if nargin<2
   error(['This function needs two arguments']);
end

fieldNames = fieldnames(dataStruct);
Nparam = size(fieldNames, 1);

%fileName = [fileName, '.txt'];
fid = fopen(fileName, 'w');
if fid == -1
  error(['Cannot open ', fileName ]);
end

fprintf(fid,'list(');
for i=1:Nparam
  fn = fieldNames(i);
  fval = fn{1};
  val = getfield(dataStruct, fval);
  [sfield1, sfield2]= size(val);

  msfield = max(sfield1, sfield2);
  newfval = strrep(fval, '_', '.');
  if ((sfield1 == 1) && (sfield2 == 1))  % if the field is a singleton
    fprintf(fid, '%s=%G',newfval, val);

  %
  % One-D array:
  %   beta = c(6, 6, ...)
  % 
  % 2-D or more:
  %   Y=structure(
  %     .Data = c(1, 2, ...), .Dim = c(30,5))
  %
  elseif ((length(size(val)) == 2) && ((sfield1 == 1) || (sfield2 == 1)))
    fprintf(fid, '%s=c(',newfval);
    for j=1:msfield
      if (isnan(val(j)))
        fprintf(fid,'NA');
      else
        % format for winbugs
        fprintf(fid,wb_strval(val(j)));
      end
      if (j<msfield)
        fprintf(fid, ', ');
      else
        fprintf(fid, ')');
      end
    end
  else
    % non-trivial 2-D or more array
    valsize    = size(val);
    alldatalen = prod(valsize);
    alldata = reshape(val', [1, alldatalen]);
    fprintf(fid, '%s=structure(.Data=c(', newfval);
    for j=1:alldatalen
      if (isnan(alldata(j)))
        fprintf(fid,'NA');
      else
        % format for winbugs
        fprintf(fid,wb_strval(alldata(j)));
      end
      if (j < alldatalen)
        fprintf(fid,',');
      else
        fprintf(fid,'), .Dim=c(', alldata(j));
      end
    end

    for j=1:length(valsize)
      if (j < length(valsize))
        fprintf(fid, '%G,', valsize(j));
      else
        fprintf(fid, '%G))', valsize(j));
      end
    end
  end
  if (i<Nparam)
    fprintf(fid, ', ');
  else
    fprintf(fid, ')\n');
  end
end
fclose(fid);

%%%%%%%%

function s = wb_strval(v)
% Converts numeric value to a string that is acceptable by winbugs.
% This is most problematic for exponent values which must have at least 1
% decimal and two digits for the exponent. So 1.0E+01 rather than 1E+001
% Note that only Matlab on PC does 3 digits for exponent.
s = sprintf('%G', v);
if strfind(s, 'E')
   if length(strfind(s, '.')) == 0
      s = strrep(s, 'E', '.0E');
   end
   s = strrep(s, 'E+0', 'E+');
   s = strrep(s, 'E-0', 'E-');
end

%%%%%%%%

function f = fullfileKPM(varargin)
% fullfileKPM Concatenate strings with file separator, then convert it to a/b/c
% function f = fullfileKPM(varargin)

f = fullfile(varargin{:});
f = strrep(f, '\', '/');

%%%%%%%%

function A = structsToArrays(S)
% Suppose S is this struct array
%
% S(c).X1(s)
% S(c).X2(s,i)
% S(c).X3(s,i,j)
%
% where s=1:N in all cases 
%
% Then we return
% A.X1(c,s)
% A.X2(c,s,i)
% A.X3(c,s,i,j)

C = length(S);
fld = fieldnames(S);
A = [];
for fi=1:length(fld)
  fname = fld{fi};
  tmp = getfield(S(1), fname);
  sz = size(tmp);
  psz = prod(sz);
  data = zeros(C, psz);
  for c=1:C
    tmp = getfield(S(c), fname);
    %data = cat(1, data, tmp);
    data(c,:) = tmp(:)';
  end
  if sz(2) > 1 % vector or matrix variable
    data = reshape(data, [C sz]);
  end
  A = setfield(A, fname, data);
end
  
 
%%%%%%%%%%%%

function [Rhat, m, s] = EPSR(samples)
% Note - not currently used.
%
% function [R, m, s] = EPSR(samples)
% "estimated potential scale reduction" statistics due to Gelman and Rubin.
% samples(i,j) for sample i, chain j
% 
% R = measure of scale reduction - value below 1.1 means converged:
%                                  see Gelman p297
% m = mean(samples)
% s = std(samples)

% This is the same as the netlab function convcalc(samples')

[n m] = size(samples);
meanPerChain = mean(samples,1); % each column of samples is a chain
meanOverall = mean(meanPerChain);

% Rhat only works if more than one chain is specified.
if m > 1
  % between sequence variace
  B = (n/(m-1))*sum( (meanPerChain-meanOverall).^2);

  % within sequence variance
  varPerChain = var(samples);
  W = (1/m)*sum(varPerChain);

  vhat = ((n-1)/n)*W + (1/n)*B;
  Rhat = sqrt(vhat/(W+eps));
else
  Rhat = nan;
end
   
m = meanOverall;
s = std(samples(:));

%%%%%%%%%

function stats = computeStats(A)

fld = fieldnames(A);
N = length(fld);
stats = struct('Rhat',[], 'mean', [], 'std', []);
for fi=1:length(fld)
  fname = fld{fi};
  samples = getfield(A, fname);
  sz = size(samples);
  clear R m s
  % samples(c, s, i,j,k)
  Nchains = sz(1);
  Nsamples = sz(2);

  st_mean_per_chain = mean(samples, 2);
  st_mean_overall   = mean(st_mean_per_chain, 1);
  
  % "estimated potential scale reduction" statistics due to Gelman and
  % Rubin.
  if Nchains > 1
    B = (Nsamples/Nchains-1) * ...
       sum((st_mean_per_chain - repmat(st_mean_overall, [Nchains,1])).^2);
    varPerChain = var(samples, 0, 2);
    W = (1/Nchains) * sum(varPerChain);
    vhat = ((Nsamples-1)/Nsamples) * W + (1/Nsamples) * B;
    Rhat = sqrt(vhat./(W+eps));
  else
    Rhat = nan;
  end

  % reshape and take standard deviation over all samples, all chains
  samp_shape = size(squeeze(st_mean_overall));
  % padarray is here http://www.mathworks.com/access/helpdesk/help/toolbox/images/padarray.html
  %reshape_target = padarray(samp_shape, [0 1], Nchains * Nsamples, 'pre');
  reshape_target = [Nchains * Nsamples, samp_shape]; % fix from Andrew Jackson  andrew@msoc.mrc.gla.ac.uk
  reshaped_samples = reshape(samples, reshape_target);
  st_std_overall = std(reshaped_samples);

  if ~isnan(Rhat)
    stats.Rhat = setfield(stats.Rhat, fname, squeeze(Rhat));
  end
  
  % special case - if mean is a 1-d array, make sure it's long
  squ_mean_overall = squeeze(st_mean_overall);
  st_mean_size = size(squ_mean_overall);
  if (length(st_mean_size) == 2) && (st_mean_size(2) == 1)
    stats.mean = setfield(stats.mean, fname, squ_mean_overall');
  else
    stats.mean = setfield(stats.mean, fname, squ_mean_overall);
  end 
  
  stats.std = setfield(stats.std, fname, squeeze(st_std_overall));
end

%%%%%%%%%%%%

function S=bugs2mat(file_ind,file_out,dir)
%BUGS2MAT  Read (Win)BUGS CODA output to matlab structure
%
% S=bugs2mat(file_ind,file_out,dir)
%  file_ind - index file (in ascii format)
%  file_out - output file (in ascii format)
%  dir      - directory where the files are found (optional)
%  S        - matlab structure, with CODA variables as fields
%
% The samples are stored in added 1'st dimension,
% so that 2 x 3 variable R with 1000 samples would be
% returned as S.R(1000,2,3)
%
% Note1: the data is returned in a structure that makes extraction
% of individual sample sequencies easy: the sequencies are
% directly Nx1 double vectors, as for example S.R(:,1,2).
% The computed statistics must, however, be squeezed,
% as mean(S.R,1) is a 1x2x2 matrix.
%
% Note2: in variable names "." is replaced with "_"

% To change the output structure, edit the 'eval' line in the m-file.
% For example, to return all samples as a cell, wich possibly varying
% number of samples for elements of a multidimensional variable,
% cange the 'eval' line to
%    eval(['S.' varname '={samples};']);
% Then the samples of R(2,1) would be returned as cell S.R(2,1)

% (c) Jouko.Lampinen@hut.fi, 2000
% 2003-01-14 Aki.Vehtari@hut.fi - Replace "." with "_" in variable names
% slightly modified by Maryam Mahdaviani, August 2005 (to suppress redundant output)

if nargin>2,
  file_ind=[dir '/' file_ind];
  file_out=[dir '/' file_out];
end

ind=readfile(file_ind);

data=load(file_out);

Nvars=size(ind,1);
S=[];
for k=1:Nvars
  [varname,indexstr]=strtok(ind(k,:));
  varname=strrep(varname,'.','_');
  indices=str2num(indexstr);
  if size(indices)~=[1 2]
     error(['Cannot read line: [' ind(k,:) ']']);
  end
  sdata = size(data);
  %indices
  samples=data(indices(1):indices(2),2);
  varname(varname=='[')='(';
  varname(varname==']')=')';
  leftparen=find(varname=='(');
  outstruct=varname;
  if ~isempty(leftparen)
     outstruct=sprintf('%s(:,%s',varname(1:leftparen-1),varname(leftparen+1:end));
  end
  eval(['S.' outstruct '=samples;']);
end

function T=readfile(filename)
f=fopen(filename,'r');
if f==-1, fclose(f); error(filename); end
i=1;
while 1
 clear line;
 line=fgetl(f);
 if ~isstr(line), break, end
 n=length(line);
 T(i,1:n)=line(1:n);
 i=i+1;
end
fclose(f);


% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_optons(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a 
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling 
% func(0, 1, 'V', 2) has the same effect.  The function call
% 
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = process_options(args, varargin)

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
  error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
  warn = 1;
  nout = n / 2;
else
  warn = 0;
  nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
  varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
  found = 0;
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
      args{i}
    else
      nunused = nunused + 1;
      unused{2 * nunused - 1} = args{i};
      unused{2 * nunused} = args{i + 1};
    end
  end
end

% Assign the unused arguments
if (~warn)
  if (nunused)
    varargout{nout} = unused;
  else
    varargout{nout} = cell(0);
  end
end
