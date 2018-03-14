function [ yt,zt ] = PLEMLambda( A,y,z,B,prior )
%PLEMLAMBDA Summary of this function goes here
%   Detailed explanation goes here
[N,K] = size(y);
[M,L] = size(z);
Lambda1 = B*(z'*z);
b1 = A*z;
temp = b1*log(Lambda1')-ones(N,L)*Lambda1';
temp = temp - repmat(mean(temp,2),[1,K]);
yt = exp(temp);
if strcmp(prior, 'prior')
    pi1 = sum(y)/N;
    yt = yt.*repmat(pi1,[N,1]);%include prior
end
yt = yt./repmat(sum(yt,2),[1,K]); %normalize

Lambda2 = B'*(y'*y);
b2 = A'*y;
temp = b2*log(Lambda2')-ones(M,K)*Lambda2';
temp = temp - repmat(mean(temp,2),[1,L]);
zt = exp(temp);
if strcmp(prior, 'prior')
    pi2 = sum(z)/M;
    zt = zt.*repmat(pi2,[M,1]);%include prior
end
zt = zt./repmat(sum(zt,2),[1,L]);
end

