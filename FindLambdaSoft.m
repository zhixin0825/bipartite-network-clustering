function [ Lambda ] = FindLambdaSoft( A,Piy,Piz )
%FINDLAMBDASOFT Summary of this function goes here
%   Detailed explanation goes here
Piy1 = bsxfun(@rdivide, Piy, sum(Piy));
%Piz1 = bsxfun(@rdivide, Piz, sum(Piz));
Lambda = Piy1'*A*Piz;
Lambda(Lambda==0) = 0.01;
end

