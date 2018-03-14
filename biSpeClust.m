function [ yh, zh, U, S, V ] = biSpeClust( A, K, L, varargin)
%SPECTRALCLUSTERING Summary of this function goes here
%   Detailed explanation goes here


%%% Optional arguments 
parser = inputParser;
addOptional(parser,'type','U')  % could be 'U' or 'US', defualt is U
addOptional(parser,'kmeans_rep',10)  
addOptional(parser,'rmax', min(K,L))  % maximum number of singular values kept
addOptional(parser,'NORMALIZE', false)  % whether to noramlize rows of the singular vector matrix
addOptional(parser,'LAPLACIAN', true) 

parse(parser, varargin{:});
type = parser.Results.type;
kmeans_rep = parser.Results.kmeans_rep;
rmax = parser.Results.rmax;
NORMALIZE = parser.Results.NORMALIZE;
LAPLACIAN = parser.Results.LAPLACIAN;
%%% end optional arumgnets

if LAPLACIAN
    %D_1 = sum(A,2);
    %D_2 = sum(A,1);
    %g1 = safe_diag_pwr(D_1);
    %g2 = safe_diag_pwr(D_2);     
    g1 = infHandle(diag(sum(A,2).^(-0.5)));
    g2 = infHandle(diag(sum(A,1).^(-0.5)));
    A_n = g1*A*g2;
else
    A_n = A;
end
   
[U,S,V] = svds(A_n,rmax);

if NORMALIZE    
    U = row_l2_normalize(U);
    V = row_l2_normalize(V);
end

% if options.verbose
%     kmopts = statset('Display','iter');
% else
%     kmopts = statset('Display','off'); 
% end
switch lower(type)
    case 'u' % default
        yh = kmeans(U, K, 'replicates', kmeans_rep, 'onlinephase','off');%,'Options',kmopts);
        zh = kmeans(V, L, 'replicates', kmeans_rep, 'onlinephase','off');%,'Options',kmopts);
    case 'us'
        yh = kmeans(U*S, K, 'replicates', kmeans_rep,'onlinephase','off');%,'Options',kmopts);
        zh = kmeans(V*S, L, 'replicates', kmeans_rep,'onlinephase','off');%,'Options',kmopts);
end

end % biSpectralClustering

function y = infHandle(x)
y = x;
y(isinf(x)) = 0;
end


function G = safe_diag_pwr(d)
    idx = d ~= 0;
    dinv = zeros(size(d));
    dinv(idx) = d(idx).^(-.5);
    G = diag(dinv);
end
