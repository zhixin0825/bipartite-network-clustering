function [ Pi ] = ytoPi( y,K )
%YTOPI Summary of this function goes here
%   Detailed explanation goes here
N = numel(y);
Pi = zeros(N,K);
for i = 1:N
    Pi(i,y(i))=1;
end

end

