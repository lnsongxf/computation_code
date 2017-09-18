clear;

% Input
input = linspace(0,1,100);
inputSize = length(input);

% Preallocate output space
output = zeros(size(input));

% Run mex
hello_mex;