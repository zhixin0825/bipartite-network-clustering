function [ yt, pi, Lambda, t ] = PLEMOneStep( A, yt, zt, varargin)
%PLEMONESTEP Summary of this function goes here
%   Detailed explanation goes here

%%% Assign optional arguments
parser = inputParser;
parser.KeepUnmatched = true;

addOptional(parser,'inner_hard',true)  
addOptional(parser,'inner_loop',true)  
addOptional(parser,'prior', true)  
addOptional(parser,'Tmax_inner', 10)  
addOptional(parser,'tol_inner', 1e-4)    % convergence tolerance
addOptional(parser,'verb', 2)  % verbose flag (levels = 0,1,2,...)

parse(parser, varargin{:});
inner_hard = parser.Results.inner_hard;
inner_loop = parser.Results.inner_loop;
prior = parser.Results.prior;
Tmax_inner = parser.Results.Tmax_inner;
verb = parser.Results.verb;
tol_inner = parser.Results.tol_inner;
%%% end optional argument

safe_exp = @(X) exp( bsxfun(@plus, X , -max(X,[],2)) );
%tau_norm = @(dtau) norm(dtau,'inf'); % this is useful for soft_labels
tau_norm = @(dtau) mean(sum((dtau).^2/2,2));

%T = 10;
[N,K] = size(yt);
L = size(zt,2);
%Y = zeros(N,T);
b = A*zt;
pi = ones(1,K)/K;

if ~inner_loop, Tmax_inner = 1; end  % if no inner loop set total numeber of iterations to 1

for t = 1:Tmax_inner
    yt_old = yt;
    
    Lambda = FindLambdaSoft(A,yt,zt);
    
    %temp = b*log(Lambda')-ones(N,L)*Lambda';
    %temp = temp - repmat(mean(temp,2),[1,K]);
    %yt = exp(temp);
    yt = safe_exp(b*log(Lambda')-ones(N,L)*Lambda');
    
    if prior %include prior
         pi = sum(yt)/N;
         yt = yt.*repmat(pi,[N,1]);
    end
    yt = row_normalize_ell1(yt);
    
    if inner_hard % tansfrom soft labels to hard
        yt = soft_to_hard(yt);
    end
    
    delta = tau_norm(yt-yt_old);
    if verb > 1, fprintf('---- inner>  %3d | del = %3.3e \n', t, delta), end
    
    if delta < tol_inner, break, end
   
end
