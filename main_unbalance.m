addpath('common')

%B = [1,2,3,4,5,6;2,3,4,5,6,1;3,4,5,6,1,2;4,5,6,1,2,3];
B = [1,2,3,4,5,6,7,8,9,10,11,12;4,5,6,7,8,9,10,11,12,1,2,3;7,8,9,10,11,12,1,2,3,4,5,6;10,11,12,1,2,3,4,5,6,7,8,9];
 
[K,L] = size(B);
%Kv = [1,1,1,1];
%Lv = [1,1,1,1,1,1];
%Kv = [1,4,6,9];
%Lv = [1,3,4,6,7,9];
Kv = ones(4,1);
Lv = ones(12,1);
beta = sqrt(max([max(Kv)/min(Kv),max(Lv)/max(Lv)])); 
%n0_vec = [50 100 200 500 1000];
%n0_len = length(n0_vec);
n0_len = 10;
n0_vec = round(logspace(log10(30), log10(2000), n0_len)); 


%%
T = 100;
acc_all = zeros(n0_len, 4, 3, T);
nmi_all = zeros(n0_len, 4, 3, T);
%nmi_all = zeros(n0_len, 4, T);
dt_all = zeros(n0_len, 4, T);

prefactor = 0.5;
al = 0.75;
%%
%for t = 1:T
parfor t = 1:T  % This is parallel for; use with caution
    dt = zeros(n0_len,4);
    %nmi = zeros(n0_len,4,2);
    nmi = zeros(n0_len,4,3);
    acc = zeros(n0_len,4,3);
    
    for j = 1:n0_len
        fprintf('--- t = %3d, j = %3d ---\n',t,j)
        n0 = n0_vec(j);
        m = L*n0;
        n = K*n0;
        compute_overll_acc = @(y,yt,z,zt) (compute_acc(y,yt)*n+compute_acc(z,zt)*m)/(n+m);
        compute_overll_nmi = @(y,yt,z,zt) compute_mutual_info(blkdiag(y,z),blkdiag(yt,zt));
        
        y = generate_random_labels(n,Kv);
        z = generate_random_labels(m,Lv);
        ny = sum(y);
        nz = sum(z);
        P = prefactor*B*(log(m*n)^al)/sqrt(m*n);
        % P = prefactor*B*log(m*n)/sqrt(m*n);
        A = genSBM2(P,y,z);
         
        opts = {'inner_loop', true, 'prior', false, 'new_y', false, ...
                'verb', 1, 'tol_inner', 1e-6, 'tol_outer', 1e-6};
        for mtd = 1:4
            tic,
            switch mtd
                case 1                
                   [ yt, zt ] = biSpeClust( A, K, L, 'LAPLACIAN', true, 'type', 'u', 'kmeans_rep',10);
                   yt = label_vec2mat(yt);
                   zt = label_vec2mat(zt);
                   yt_init = yt;
                   zt_init = zt;
                case 2
                    [yt, zt] = PLEMLambda( A,y,z,P,'noprior' );
                    yt = soft_to_hard(yt);
                    zt = soft_to_hard(zt);
                    %PLLLH2 = PLLH(A, yt, zt);
                case 3
                    [yt, zt] = PLEM( A, yt_init, zt_init,  'inner_hard', true, 'outer_hard', true, opts{:});
                case 4
                    [yt, zt] = PLEM( A, yt_init, zt_init, 'inner_hard', false, 'outer_hard', false, opts{:});
                    
            end       
            dt(j,mtd) = toc;
        
            nmi(j,mtd,1) = compute_mutual_info(yt,y);
            nmi(j,mtd,2) = compute_mutual_info(zt,z);
            nmi(j,mtd,3) = compute_overll_nmi(y, yt, z, zt);
            acc(j,mtd,1) = compute_acc(y, yt);
            acc(j,mtd,2) = compute_acc(z, zt);
            acc(j,mtd,3) = compute_overll_acc(y, yt, z, zt);
        end
        
%         tic,
%         [yt2,zt2] = PLEMLambda( A,y,z,P,'noprior' );
%         yt2 = soft_to_hard(yt2);
%         zt2 = soft_to_hard(zt2);
%         nmi(j,2,1) = compute_mutual_info(yt2,y);
%         nmi(j,2,2) = compute_mutual_info(zt2,z);
%         nmi(j,2,3) = compute_overll_nmi(y, yt2, z, zt2);
%         acc(j,2) = compute_overll_acc(y, yt2, z, zt2);
%         PLLLH2 = PLLH(A,yt2,zt2);
%         [yt3,zt3] = PLEM( A,yt1,zt1, 'innerhard', 'outerhard', 'noinnerloop', 'noprior', 'oldy' );
%         dt(j,2) = toc;
     
        
        
%         tic
%         [yt3, zt3] = PLEM( A, yt1, zt1,  'inner_hard', true, 'outer_hard', true, opts{:});
%         nmi(j,3,1) = compute_mutual_info(yt3,y);
%         nmi(j,3,2) = compute_mutual_info(zt3,z);
%         nmi(j,3,3) = compute_overll_nmi(y, yt3, z, zt3);
%         acc(j,3) = compute_overll_acc(y, yt3, z, zt3);
%         dt(j,3) = toc;
       
%         tic
%         [yt4, zt4] = PLEM( A, yt1, zt1, 'inner_hard', false, 'outer_hard', false, opts{:});
%         [yt4,zt4] = PLEM( A,yt1,zt1, 'innersoft', 'outersoft', 'noinnerloop', 'noprior', 'oldy' );
%         yt4 = soft_to_hard(yt4);
%         zt4 = soft_to_hard(zt4);
%         nmi(j,4,1) = compute_mutual_info(yt4,y);
%         nmi(j,4,2) = compute_mutual_info(zt4,z);
%         nmi(j,4,3) = compute_overll_nmi(y, yt4, z, zt4);
%         acc(j,4) = compute_overll_acc(y, yt4, z, zt4);
%         dt(j,4) = toc;
        
      
    end

    nmi_all(:,:,:,t) = nmi;
    acc_all(:,:,:,t) = acc;
    dt_all(:,:,t) = dt;
%     java.io.File.createTempFile(sprintf('%3d',t), 'temp');
end

%%
%fprintf('\nrunning times = %s \n', sprintf('%4.2f  ', mean(dt_all,2)) )

%%
nmi_avg = mean(nmi_all,4);
acc_avg = mean(acc_all,4);
dt_avg = mean(dt_all,3);
result_fname = strrep(sprintf('results_C%2.2f_a%2.2f_b%2.2f_T%d_K%d_L%d',prefactor, al,beta, T,K,L),'.','p');
save(sprintf('%s.mat',result_fname))

%%
%load('results.mat')
title_str = sprintf('C = %2.2f, \\alpha = %2.2f, \\beta=%2.2f, K=%d, L=%d',prefactor,al,beta,K,L);
%title_str = sprintf('C = %2.2f, \\alpha = %2.2f',prefactor,al);
figure(1), clf,
%subplot(131),hold on
%colors = parula(4);
colors = get(gca,'ColorOrder');
markers = {'-.',':','--','-'};
for i = 1:4
    %h(i) = plot_ci_bands(n0_vec, squeeze(nmi_all(:,i,3,:)), colors(i,:), @plot);
    %h(i) = semilogx(n0_vec, nmi_avg(:,i,1),'LineWidth',2,'color',colors(i,:)), hold on
    h(i) = plot(n0_vec, nmi_avg(:,i,1),  markers{i}, ...
        'LineWidth',2,'color',colors(i,:)); hold on
end

lgd = legend(h,{'Spectral','Oracle','Hard', 'Soft'}, ...
    'Position',[0.65 .15 0.2 0.2]);
legend('boxoff')
%axis([0,N+M,0,1]);
xlabel('# of nodes per cluster')
ylabel('NMI (Overall)')
title(title_str,'FontWeight','Normal')
axis tight
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4.5];
fig.PaperPositionMode = 'manual';
print('-depsc2',sprintf('%s_nmi.eps',result_fname))

%%
figure(2), clf, 
for i = 1:4  
    %plot_ci_bands(n0_vec, squeeze(acc_all(:,i,3,:)), colors(i,:), @plot)    
    semilogx(n0_vec, log(1-acc_avg(:,i,3)),markers{i},'LineWidth',2,'color',colors(i,:)), hold on
end
ylabel(' log miss. (overall)')
xlabel('# of nodes per cluster')
axis tight
title(title_str,'FontWeight','Normal')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 5 4.5];
fig.PaperPositionMode = 'manual';
print('-depsc2',sprintf('%s_miss.eps',result_fname))
