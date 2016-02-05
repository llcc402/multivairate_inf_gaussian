% Input:
%      data        a matrix. The rows are corresponding to the observations 
%                  and the columns are the variables.
%      alpha       a scalar. The concentration parameter of the DP.
%      sigma       a scalar. The precision of the prior of the mean
%                  of the gaussians.
%      gamma       a scalar. The shape and the rate parameter of a Gamma 
%                  distribution, which is the prior of the precision of the
%                  Gaussians.
%      actN        a scalar. The maximum number of the activated positions
%                  of the DP.
%      maxIter     a scalar. The maximum number of the Gibbs iterations.
% Output:
%      z           a row vector. The length of z is the number of the
%                  observations and the values are the indicators of the
%                  observations.
%      centers     a matrix. The rows are the centers and the columns are
%                  the variables.
function [z,centers, G0, precisions] = inf_gaussian(data, alpha, sigma, gamma, actN, maxIter)
if nargin < 5
    actN = 100;
end
if nargin < 6
    maxIter = 1000;
end

%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
% the indicators
z = ones(1, size(data, 1));

% parameters of the Gaussians
centers = zeros(actN, size(data, 2));
precisions = zeros(actN, size(data, 2));

% the weights of the DP
G0 = zeros(1, actN);

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    
    % aggregate the indicators
    ix = accumarray(z', 1:length(z), [], @(x){x});
    
    % number of observations taking on specific atom
    num_atom = zeros(1, length(ix));
    
    % sample positions of G0
    for i = 1:length(ix)
        if isempty(ix{i})
            num_atom(i) = 0;
            % sample a gaussian from the prior
            precisions(i,:) = gamrnd(gamma * ones(1,size(data, 2)), ...
                                     ones(1,size(data, 2)) / gamma ...
                                     );
            centers(i,:) = mvnrnd(zeros(1, size(data, 2)), ...
                                  eye(size(data, 2)) / sigma ...
                                  );
        else
            num_atom(i) = length(ix{i});
            % sample a gaussian from the posterior
            a = gamma + 1/2;
            exp_error = sum((data(ix{i},:) - repmat(centers(i,:), length(ix{i}), 1)) .^ 2);
            b = gamma + exp_error / 2;
            precisions(i,:) = gamrnd(a, b .^ (-1));
            
            centers(i,:) = sum(data(ix{i},:)) .* precisions(i,:) ./ (length(ix{i}) * precisions(i,:) + sigma);
        end
    end
    
    if length(ix) < actN
        precisions(length(ix)+1:end, :) = gamrnd(gamma * ones(actN - length(ix), size(data,2)),...
                                         ones(actN - length(ix), size(data, 2)) / gamma...
                                         );
        centers(length(ix)+1:end, :) = mvnrnd(zeros(1, size(data, 2)), ...
                                      eye(size(data, 2)) / sigma, ...
                                      actN - length(ix) ...
                                      );
    end
    
    % sample weights of G0
    if length(num_atom) < actN
        counts = [num_atom, zeros(1, actN - length(num_atom))];
    else
        counts = num_atom;
    end
    a = counts + 1;
    b = [cumsum(counts(2:end), 'reverse'), 0];
    b = b + alpha;
    V = betarnd(a,b);
    G0 = V;
    V = cumprod(1-V);
    G0(2:end) = G0(2:end) .* V(1:end-1);
    
    % sample z
    first = sum(log(precisions), 2);
    for i = 1:size(data, 1) 
        second = sum((repmat(data(i,:), actN, 1) - centers) .^ 2 .* precisions, 2);
        log_likelihood = first - second;
        log_post = log_likelihood' + log(G0);
        log_post = log_post - max(log_post);
        post = exp(log_post);
        post = post / sum(post);
        [~, ~, z(i)] = histcounts(rand(1), [0, cumsum(post)]);
    end
end

end
    














