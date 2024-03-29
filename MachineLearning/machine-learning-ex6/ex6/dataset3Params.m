function [C, sigma] = dataset3Params(X, y, Xval, yval)
%DATASET3PARAMS returns your choice of C and sigma for Part 3 of the exercise
%where you select the optimal (C, sigma) learning parameters to use for SVM
%with RBF kernel
%   [C, sigma] = DATASET3PARAMS(X, y, Xval, yval) returns your choice of C and 
%   sigma. You should complete this function to return the optimal C and 
%   sigma based on a cross-validation set.
%

% You need to return the following variables correctly.
C = 1;
sigma = 0.3;

% ====================== YOUR CODE HERE ======================
% Instructions: Fill in this function to return the optimal C and sigma
%               learning parameters found using the cross validation set.
%               You can use svmPredict to predict the labels on the cross
%               validation set. For example, 
%                   predictions = svmPredict(model, Xval);
%               will return the predictions on the cross validation set.
%
%  Note: You can compute the prediction error using 
%        mean(double(predictions ~= yval))
%
C_vec = [0.01 0.03 0.1 0.3 1 3 10 30];
sigma_vec = [0.01 0.03 0.1 0.3 1 3 10 30];
predictions_error = zeros(length(C_vec),length(sigma_vec));
precision = zeros(length(C_vec),length(sigma_vec));
recall = zeros(length(C_vec),length(sigma_vec));
F1_score = zeros(length(C_vec),length(sigma_vec));

for i = [1:length(C_vec)]
  for j = [1:length(sigma_vec)]
    C = C_vec(i);
    sigma = sigma_vec(j);
    kernelFunction = @(x1,x2)gaussianKernel(x1,x2,sigma);
    model = svmTrain(X,y,C,kernelFunction);
    predictions = svmPredict(model,Xval);
    predictions_error(i,j) = mean(double(predictions!=yval));
    precision(i,j) = sum(predictions==yval)/sum(predictions);
    recall(i,j) = sum(predictions==yval)/sum(yval);
    F1_score(i,j) = 2*precision(i,j)*recall(i,j)/(precision(i,j)+recall(i,j));
  endfor
endfor

%F1_score = 2*precision.*recall./(precision+recall);
precision
recall
F1_score
predictions_error

[M I] = min(predictions_error(:));
[I_row,I_col] = ind2sub(size(predictions_error),I);
C = C_vec(I_row);
sigma = sigma_vec(I_col);

% =========================================================================

end
