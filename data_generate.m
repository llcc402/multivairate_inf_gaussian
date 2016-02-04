function data = data_generate(n)
if nargin < 1
    n = 1000;
end

theta = rand(n, 1) * 2 * pi;
r = 1;
epsilon = .2;
x = r * cos(theta) + rand(n, 1) * epsilon;
y = r * sin(theta) + rand(n, 1) * epsilon;

data1 = [x, y];

theta = rand(n, 1) * 2 * pi;
r = 2;
x = r * cos(theta) + rand(n, 1) * epsilon;
y = r * sin(theta) + rand(n, 1) * epsilon;

data2 = [x, y];

data = [data1; data2];
plot(data(:,1), data(:,2), 'o')