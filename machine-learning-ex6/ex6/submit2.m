x1 = sin(1:10)';
x2 = cos(1:10)';
ec = 'the quick brown fox jumped over the lazy dog';
wi = 1 + abs(round(x1 * 1863));
wi = [wi ; wi];

sim = gaussianKernel(x1, x2, 2);
out = sprintf('%0.5f ', sim);

load('ex6data3.mat');
    [C, sigma] = dataset3Params(X, y, Xval, yval);
    out = sprintf('%0.5f ', C);
    out = [out sprintf('%0.5f ', sigma)];