addpath('common')
B = [1,2,3,4,5,6;
     2,3,4,5,6,1;
     3,4,5,6,1,2;
     4,5,6,1,2,3];
[K,L] = size(B);



n0 = 500;
m = L*n0;
n = K*n0;
y = generate_random_labels(n,K);
z = generate_random_labels(m,L);

ny = sum(y);
nz = sum(z);
P = .5*B*log(m*n)/sqrt(m*n);
A = genSBM2(P,y,z);
fprintf('density = %3.3f\n', nnz(A)/(n*m))

figure(2), clf, spy(A)

M = y*P*z';

[ ys, zs, Us ] = biSpeClust( M, K, L, 'LAPLACIAN', true, 'type', 'u');%, 'rmax',5);
[ yh, zh, U ] = biSpeClust( A, K, L, 'LAPLACIAN', true, 'type', 'u');%, 'kmeans_rep',50, 'rmax', 3);
%[ yh1, U1 ] = SpectralClustering(A,min(K,L),K,'U');

%%
figure(1), clf
% idx = 3:5;
% U = U(:,idx);
subplot(131), scatter3(Us(:,1),Us(:,2),Us(:,3),'.')
axis tight
subplot(132), scatter3(U(:,1),U(:,2),U(:,3),'.')
axis tight
%U1 = U1(:,idx);
colors = {'r','g','b','y'};

subplot(133), hold on
%yh_lab = kmeans(U(:,1:3),4); 
yh_lab = label_mat2vec(yh);
for k = 1:4
    idx = yh_lab == k;
    scatter3(U(idx,1),U(idx,2),U(idx,3),[colors{k} '.'])
    grid on
    axis tight
    view(3)
    
end

fprintf('%f, %f\n', compute_mutual_info(ys,y),compute_mutual_info(zs,z))
fprintf('%f, %f\n', compute_mutual_info(yh,y),compute_mutual_info(zh,z))
fprintf('%f\n', compute_mutual_info(yh1,y))

fprintf('acc = %f\n', compute_acc(yh1,y))
%%
% yt1 = SpectralClustering(A,min(K,L),K,'U');
%compute_mutual_info(yt1,yh)
