function [output] = sumExp(vars,data)
%SUMEXP Summary of this function goes here
%   Detailed explanation goes here
output = 0;
for cntI = 1:length(vars)/2
     output = output + vars(cntI*2-1).*exp(-(data)./vars(cntI*2));
end
end

