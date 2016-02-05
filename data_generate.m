% Input:
%      n      a scalar. The number of observations in each cluster.
%      k      a scalar. The number of clusters.
%      d      a scalar. The dimension of the data.
% Output:
%      data   a matrix of order (n * k) times d. 
function [data, centers] = data_generate(n, k, d)
if nargin < 1
    n = 500;
end
if nargin < 2
    k = 5;
end
if nargin < 3
    d = 2;
end

%--------------------------------------------------------------------------
% STEP 1: Generate centers
%--------------------------------------------------------------------------
sigma = 40;
centers = mvnrnd(zeros(1, d), eye(d) * sigma, k);

%--------------------------------------------------------------------------
% STEP 2: Generate mixing measure
%--------------------------------------------------------------------------
mixing = rand(1, k);
mixing = mixing / sum(mixing);

%--------------------------------------------------------------------------
% STEP 3: Generate latent variable
%--------------------------------------------------------------------------
r = rand(1, n*k);
[~, ~, z] = histcounts(r, [0, cumsum(mixing)]);

%--------------------------------------------------------------------------
% STEP 3: Generate data 
%--------------------------------------------------------------------------
data = zeros(n*k, d);
for i = 1 : n*k
    data(i,:) = mvnrnd(centers(z(i),:), eye(d));
end

if d == 2
    plot(data(:,1), data(:,2), '.')
end

end