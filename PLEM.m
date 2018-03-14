function [ yt, zt, t ] = PLEM( A, yt, zt, varargin )
%function [ y,z,t ] = PLEM( A,yt,zt, innerhard, outerhard, innerloop, prior, newy )
%PLEM Summary of this function goes here
%   Detailed explanation goes here

%%% Assign optional arguments
parser = inputParser;
parser.KeepUnmatched = true;

addOptional(parser,'outer_hard', true)  
addOptional(parser,'new_y',  true)    
addOptional(parser,'Tmax_outer', 100)  
addOptional(parser,'tol_outer', 1e-4)    % convergence tolerance
addOptional(parser,'verb', 2)  % verbose flag (levels = 0,1,2,...)

parse(parser, varargin{:});
outer_hard = parser.Results.outer_hard;
new_y = parser.Results.new_y;
Tmax_outer = parser.Results.Tmax_outer;
verb = parser.Results.verb;
tol_outer = parser.Results.tol_outer;
%%% end optional argument

%tau_norm = @(dtau) norm(dtau,'inf'); % this is useful for soft_labels
tau_norm = @(dtau) mean(sum((dtau).^2/2,2));

delta = zeros(1,2);
for t = 1:Tmax_outer
    yt_old = yt;
    zt_old = zt;
    
    yt = PLEMOneStep(A, yt, zt, varargin{:});
    
    if new_y
        zt = PLEMOneStep(A', zt, yt, varargin{:});
    else
        zt = PLEMOneStep(A', zt, yt_old, varargin{:});
    end
    
    if outer_hard % tansfrom soft labels to hard
        yt = soft_to_hard(yt);
        zt = soft_to_hard(zt);
    end
    
    
    delta(1) = tau_norm(yt-yt_old);
    delta(2) = tau_norm(zt-zt_old);
    if verb > 0 && mod(t,10) == 0, fprintf('outer>  %3d | del = %3.3e, %3.3e\n', t, delta(1), delta(2)), end
    
    if max(delta) < tol_outer, break, end
end
