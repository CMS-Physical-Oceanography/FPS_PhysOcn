function [output,flag] = optics_filter(input,perc)
%[output,flag] = optics_filter(input)
%
% A very rough first difference filter that just removes "spikes" and
% preserves 99% of the data (perc = 0.99)
%
todo = input;

%Trying out a simple first difference filter (This part of the code is from percentileAS.m):
X = abs(diff(todo));
% Retain 99% of values
test3 = sort(X);
val = test3(floor(perc*length(test3)));

thr = val;
aa = find(abs(diff(todo))>=thr);

output = input;
output(aa) = NaN;

flag = isnan(output);