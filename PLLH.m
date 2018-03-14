function [ llh ] = PLLH( A,y,z )
%PLLH Summary of this function goes here
%   Detailed explanation goes here

llh = halfpllh(A,y,z)+halfpllh(A',z,y);

function [ llh ] = halfpllh(A,y,z)
    Lambda = FindLambdaSoft(A,y,z);
    b = A*z;
    c = y*Lambda;
    bv = reshape(b,[],1);
    cv = reshape(c,[],1);
    llh = sum(cv.*log(bv)-cv-factorial(bv));
end
end

