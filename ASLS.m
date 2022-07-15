%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this code is based on Eilers, Boelens 'Baseline Correction with 
% Asymmetric Least Squares Smoothing' 2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function baseline_est = ASLS(data, ASLS_param)
% Inputs:
% data = nx1 column vector of data values (any units, any sample rate)
% ASLS_param >> lambda (default 1e5), p (default 0.01), noise_margin (default 0), max_iter (default 5)
%	if the struct doesn't include all fields, it will use default values for those params
%	if true, uses default parameters for baseline subtraction
%	defaults to true
% 
% Outputs:
% baseline_est = nx1 column vector of the fitted model for the baseline (same units & sample rate as data)
% 
% notes on ASLS parmaters:
% * anything above the noise margin gets weighted by p
% * anything within the noise margin gets weighted by 1-p
%     >> so 0>p>1
%     >> you want p as low as possible, but you need the fit to converge
% * lambda is a smoothing parameter that enforces the smoothness of the baseline
%     >> larger lambda ==> more smooth / less curvature (approaching flat)
%     >> smaller lambda ==> more wiggly-ness allowed in the baseline
%     -  this is maybe not unit-independent
%     -  it might approx. scale w/ how wide the events are?
%     -  if your drift is very slow, you can have a very high lambda
% * over each iteration, the algorithm is down-weighting data points outside the noise margin from the previous fit
%     >> check if converged by increasing max_iter
%     >> if use a smaller p, you might need to increase max_iter

%% Parse inputs

lambda = 1e5; % default lambda
p = 0.01; % default p
max_iter = 5; % default max_iter
noise_margin = 0; % default noise margin
if nargin == 1 || isempty(ASLS_param) || (islogical(ASLS_param) &&  ASLS_param == true)
    ASLS_param = struct('lambda', lambda, 'p', p, 'max_iter', max_iter, 'noise_margin', noise_margin); % use defaults
elseif nargin == 2
    if ~isfield(ASLS_param, 'lambda') || isempty(ASLS_param.lambda)
        ASLS_param.lambda = lambda;
    end
    if ~isfield(ASLS_param, 'p') || isempty(ASLS_param.p)
        ASLS_param.p = p;
    end
    if ~isfield(ASLS_param, 'max_iter') || isempty(ASLS_param.max_iter)
        ASLS_param.max_iter = max_iter;
    end
    if ~isfield(ASLS_param, 'noise_margin') || isempty(ASLS_param.noise_margin)
        ASLS_param.noise_margin = noise_margin;
    end    
end

%% Perform algorithm

lambda = ASLS_param.lambda;
p = ASLS_param.p;
max_iter = ASLS_param.max_iter;
noise_margin = ASLS_param.noise_margin;

N = length(data);
D = diff(speye(N), 2);
w = ones(N, 1);

for i = 1:max_iter
    W = spdiags(w, 0, N, N);
    C = chol(W + lambda * (D.' * D));
    baseline_est = C \ (C.' \ (w .* data));    
    w = p .* (data > (baseline_est+noise_margin)) + (1 - p) .* (data < (baseline_est+noise_margin));
end

end