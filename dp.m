function [z, G0, centers] = dp(data, alpha, lambda, actN, maxIter)
if nargin < 4
    actN = 100;
end
if nargin < 5
    maxIter = 1000;
end
%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
G0 = mvnrnd(zeros(1, size(data, 2)), eye(size(data, 2)), actN);
z = ones(1, size(data, 1));
centers = zeros(actN, size(data, 2));

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    
    % sample G0
    counts = histcounts(z, 1:actN+1);
    a = counts + 1;
    b = [cumsum(counts(2:end), 'reverse'), 0];
    b = b + alpha;
    V = betarnd(a, b);
    G0 = V;
    V = cumprod(1 - V);
    G0(2:end) = G0(2:end) .* V(1:end-1);
    
    % sample centers
    ix = accumarray(z', 1:length(z), [], @(x){x});
    for i = 1:length(ix)
        if isempty(ix{i})
            centers(i,:) = mvnrnd(zeros(1, size(data, 2)), eye(size(data, 2)) / lambda);
        else
            centers(i,:) = mean(data(ix{i}, :));
        end
    end
    if length(ix) < actN
        for i = length(ix)+1:actN
            centers(i,:) = mvnrnd(zeros(1, size(data, 2)), eye(size(data, 2)) / lambda);
        end
    end
    
    % sample z
    for i = 1:length(z)
        log_likelihood = repmat(data(i,:), actN, 1) - centers;
        log_likelihood = - sum(log_likelihood .^2, 2);
        log_post = log(G0) + log_likelihood';
        log_post = log_post - max(log_post);
        post = exp(log_post);
        post = post / sum(post);
        [~,~, z(i)] = histcounts(rand(1), [0, cumsum(post)]);
    end
end

end
    











