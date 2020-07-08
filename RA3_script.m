function [] = RA3_script(data_name, bulk_name, peak_selection_k, K2, K3, multi_donor)

rng(23);
% clc
% clear


if multi_donor == true
	task_name = sprintf('Our_%s_peak%03d_dim%d_multidonor',data_name,100*peak_selection_k,K2);
else
	task_name = sprintf('Our_%s_peak%03d_dim%d',data_name,100*peak_selection_k,K2);
end
res_name = sprintf('../resultForComparison/%s.mat',task_name);

fprintf('%s\n',task_name);
if (exist(res_name,'file') )
	fprintf('File exist\n');
	return
end

res_path = sprintf('../resultForComparison/results/%s_%s',task_name,datestr(now,30));
mkdir(res_path);

%% Load data
load(sprintf('../../supplementary_code/data/dataForComparison/%s_bulk_mat.mat',bulk_name))
load(sprintf('../../supplementary_code/data/dataForComparison/%s.mat',data_name))
if size(bulk_mat,1) > size(bulk_mat,2) % to check the Dim, may be removed in future version
	bulk_mat = bulk_mat';
end
bulk_mat = double(bulk_mat);
Y = double(count_mat); clear count_mat
cell_label = string(label_mat); clear label_mat

if multi_donor == true
	load(sprintf('../../supplementary_code/data/dataForComparison/%s_donor_labels.mat',data_name)) % labels of donors
else
	donor_label = cell_label;
end
if size(donor_label,1) < size(donor_label,2) % to check the Dim, may be removed in future version
	donor_label = donor_label';
end


%% filter features
non_zero_peak_ind = sum(Y') > 0;
Y = Y(non_zero_peak_ind,:);
bulk_mat = bulk_mat(:,non_zero_peak_ind);

filter_peak = sum(Y'>=1) >= round(size(Y,2)*0.01*peak_selection_k);
Y = Y(filter_peak,:);
bulk_mat = bulk_mat(:,filter_peak);

%% Normalize the data with TF-IDF
nfreqs = Y ./ repmat(sum(Y,1),[size(Y,1) 1]);
Y_mat = nfreqs .* repmat(log(1 + size(Y,2) ./ sum(Y,2)), [1 size(Y,2)]);

%% Visualization of Real Data
[coeff,score,~,~,explained] = pca(bulk_mat,'Centered',true,'economy',true);
bulk_components = 30;
if size(bulk_mat,1) > bulk_components
    coeff = coeff(:,1:bulk_components);
    score = score(:,1:bulk_components);
end

figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.7]);
subplot(1,2,1);
bulk_tsne = tsne((coeff' * Y_mat)' );
% size(bulk_tsne)
gscatter(bulk_tsne(:,1),bulk_tsne(:,2),cell_label);
title('tsne Bulk');

[U1, S1, V1] = fsvd(Y_mat', 10, 2);
norm_V1 = sqrt(diag(V1'*V1))';
coeff_pca = V1 ./ norm_V1;
score_pca = U1 * S1 .* norm_V1; clear U1 V1 S1;
% [coeff_pca, score_pca,~,~,explained_pca] = pca(Y_mat','Centered',true,'economy',true);
subplot(1,2,2);
% size(score_pca)
pca_tsne = tsne(score_pca(:,1:10));
gscatter(pca_tsne(:,1),pca_tsne(:,2),cell_label);
title('pca');
saveas(gcf,sprintf('%s/bulk_pca_tsne.png',res_path));


%% parameters setting

tic;

[p, n] = size(Y);
K1 = size(coeff,2);

center_Y = Y_mat; % Y_mat is the scaled data, Y is the raw count data
latent_h = coeff' * center_Y;% k*n

% standardiz
V_beta   = sqrt(var(latent_h')); %
init_V   = repmat(V_beta, p, 1) .* coeff; % p*k

%% good start of the K2 components
residual = center_Y - coeff * (coeff' * center_Y); 
[U1_res, S1_res, V1_res] = fsvd(residual', K2, 2);
norm_V1_res = sqrt(diag(V1_res'*V1_res))';
coeff_res = V1_res ./ norm_V1_res;
score_res = U1_res * S1_res .* norm_V1_res; clear U1_res V1_res S1_r
% [coeff_res,score_res,~,~,explained_res] = pca(residual','Centered',true,'economy',true);

score_rotate1 = rotatefactors(score_res(:,1:K2),'Method','varimax','Maxit',5000);
RM = linsolve(score_res(:,1:K2), score_rotate1);
coeff_res10 = coeff_res(:,1:K2);
coeff_res10_RM = coeff_res10*RM;
residual_coeff_res10_RM = residual'*coeff_res10_RM;

H2_ini = residual_coeff_res10_RM*diag(1./sqrt(var(residual_coeff_res10_RM))); % n by K
W2_ini = coeff_res10_RM*diag(sqrt(var(residual_coeff_res10_RM)));

%% Sigma_setting
center_Y_stand = center_Y - repmat(mean(center_Y,2),[1,n]);
residual_stand = center_Y_stand - coeff * (coeff' * center_Y_stand); 

[U1_res_stand, S1_res_stand, V1_res_stand] = fsvd(residual_stand', 20, 2);
norm_V1_res_stand = sqrt(diag(V1_res_stand'*V1_res_stand))';
coeff_res_stand = V1_res_stand ./ norm_V1_res_stand;
score_res_stand = U1_res_stand * S1_res_stand .* norm_V1_res_stand; 
clear U1_res_stand V1_res_stand S1_res_stand;
% [coeff_res_stand,score_res_stand,~,~,explained_res_stand] = pca(residual_stand','Centered',true,'economy',true);

score_rotate2 = rotatefactors(score_res_stand(:,1:20),'Method','varimax','Maxit',5000);
RM_stand = linsolve(score_res_stand(:,1:20), score_rotate2);
coeff_res10_stand = coeff_res_stand(:,1:20);
coeff_res10_RM_stand = coeff_res10_stand*RM_stand;
residual_coeff_res10_RM_stand = residual_stand'*coeff_res10_RM_stand;

% figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.7]);
% for i = 1:5
%     subplot(2,3,i);
%     % 	scatter(1:n,residual_coeff_res10_RM(:,i),100,cell_color,'.');
%     gscatter(1:n,score_rotate2(:,i),cell_label);
%     legend_s = findobj('type','legend');
%     delete(legend_s)
%     title(i);
% end
% saveas(gcf,sprintf('%s/standardized_varimax15.png',res_path));
% figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.7]);
% for i = 1:5
%     subplot(2,3,i);
%     % 	scatter(1:n,residual_coeff_res10_RM(:,i),100,cell_color,'.');
%     gscatter(1:n,score_rotate2(:,i+5),cell_label);
%     legend_s = findobj('type','legend');
%     delete(legend_s)
%     title(i);
% end
% saveas(gcf,sprintf('%s/standardized_varimax610.png',res_path));
% 
% 
% figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.7]);
% for i = 1:5
%     subplot(2,3,i);
%     % 	scatter(1:n,residual_coeff_res10_RM(:,i),100,cell_color,'.');
%     gscatter(1:n,score_rotate2(:,10+i),cell_label);
%     legend_s = findobj('type','legend');
%     delete(legend_s)
%     title(i);
% end
% saveas(gcf,sprintf('%s/standardized_varimax1115.png',res_path));
% 
% figure('Visible','off','units','normalized','position',[0.1,0.1,0.9,0.7]);
% for i = 1:5
%     subplot(2,3,i);
%     % 	scatter(1:n,residual_coeff_res10_RM(:,i),100,cell_color,'.');
%     gscatter(1:n,score_rotate2(:,15+i),cell_label);
%     legend_s = findobj('type','legend');
%     delete(legend_s)
%     title(i);
% end
% saveas(gcf,sprintf('%s/standardized_varimax1520.png',res_path));

H2_ini_stand = residual_coeff_res10_RM_stand*diag(1./sqrt(var(residual_coeff_res10_RM_stand))); % n by K
W2_ini_stand = coeff_res10_RM_stand*diag(sqrt(var(residual_coeff_res10_RM_stand)));

%% Select most deviant components index as varimax_sidx
varimax_sidx = [];

varimax_iidx = setdiff(1:20,varimax_sidx);
W2_sig = W2_ini_stand(:,varimax_sidx);
H2_sig = H2_ini_stand(:,varimax_sidx)';
W3_sig = W2_ini_stand(:,varimax_iidx);
H3_sig = H2_ini_stand(:,varimax_iidx)';
epsilon_hat = residual_stand - W2_sig * H2_sig - W3_sig * H3_sig;
epsilon_stard = epsilon_hat - repmat(mean(epsilon_hat,2),[1,n]);
sigma_setting = trace(epsilon_stard' * epsilon_stard)/(n*p);


%% ARD      
% Parameter Setting
sigma = 1;
sigma1 = 0.9;
sigma2 = 5;
theta = 0.1;
K = K1 + K2 + K3;
para_num = p * K;

set_sparse=true;


if multi_donor == true
    q = size(unique(donor_label),1);
    X = zeros(q, n);% 1*n
    M = containers.Map(unique(donor_label), 1:q);
    labels = cell2mat(values(M,donor_label));
    for i = 1:n
        X(labels(i),i) = 1;
    end
    Beta = zeros(p,q);
else
    q = 1;X = ones(1,n);Beta = zeros(p,q); % intercept
end


% Initialization
W1_ini = init_V(:,1:K1);
H1_ini = latent_h(1:K1,:)'*diag(1./V_beta(1:K1));
varmax_sort = [1:K2];
varmax_ind = varmax_sort(1:(K2));

H_PCA = [H1_ini H2_ini(:,varmax_ind) randn(n, K3)]';
W_PCA = [W1_ini W2_ini(:,varmax_ind) randn(p,K3)];

A_PCA = diag([ones(1,K2) ones(1,K3)]);
fix_A = false;

Gamma_PCA = ones(K,n);
Gamma_PCA((K1+1):(K2+K1),:) = zeros(K2,n);
V_PCA = [init_V(:,1:K1) zeros(p,K2+K3)];
init_time = double(toc);
% Estimation
[W_train,H_train,Gamma_train,lgp_train,alpha_train,Beta_train,sigma_train] = RA3_func(donor_label,cell_label,Y_mat,theta,sigma,sigma1,sigma2,K1,K2,K3,Gamma_PCA,A_PCA,W_PCA,V_PCA,sigma_setting,H_PCA,X,Beta,res_path,set_sparse,fix_A);
% fprintf('\n[K2 = %d, K3 = %d]: run finished. PLOTTING...\n', K2, K3);
[indicator_o P_o] = ttest(H_train(K1+1:K1+K2,:)');
if sum(indicator_o==1) == 0
    sparse_index_left = 0;
    K2_left = 0;
    H2_trun = 0;
else
    sparse_index_left = find(indicator_o == 1);
    K2_left = length(sparse_index_left);
    H2_trun = H_train(K1+sparse_index_left, :);
    for k = 1:K2_left
        TEP = H2_trun(k,:);
        TEP(TEP >= quantile(TEP,0.95)) = quantile(TEP,0.95);
        TEP(TEP <= quantile(TEP,0.05)) = quantile(TEP,0.05);
        H2_trun(k,:) = TEP;
    end
end

total_time = double(toc);

%% save results for plotting 
sc_pca_10 = score_pca(:,1:10);
bulk_proj = (coeff' * Y_mat)';
cell_colors = cell_label;
cell_labels = cell_label;
donor_colors = donor_label;
donor_labels = donor_label;
W_hat = W_train;
H_hat = H_train;
Gamma_hat = Gamma_train;
lgp = lgp_train;
A_hat = alpha_train;
Beta_hat = Beta_train;
sigma_s_o = sigma_train;
save(res_name,...
    'sc_pca_10', 'bulk_proj','cell_colors','cell_labels','donor_colors','donor_labels',...
    'W_hat','H_hat','Gamma_hat','A_hat','Beta_hat','sigma_s_o','lgp','K1','K2','K3',...
    'p', 'n', 'res_path','init_time','total_time','sparse_index_left','K2_left','H2_trun');

%% Check W, H and A estimate
alpha_H2_train = diag(alpha_train(1 : K2 ,1 : K2 ));
alpha_H3_train = diag(alpha_train((K2+1) : end,(K2+1) : end));

scale_W = sqrt(diag(W_train' * W_train))'; 
W_train_stand = W_train ./ scale_W;
H_train_stand = H_train .* scale_W';
H_train_norm = diag(H_train_stand * H_train_stand');
W_train_norm = diag(W_train_stand' * W_train_stand);

W1_train_norm = W_train_norm(1:K1);
W2_train_norm = W_train_norm((1+K1) : (K1+K2));
W3_train_norm = W_train_norm((1+K1+K2) : K);

H1_train_norm = H_train_norm(1:K1);
H2_train_norm = H_train_norm((1+K1) : (K1+K2));
H3_train_norm = H_train_norm((1+K1+K2) : K);


H_hat5 = H_train'; % n by K

figure('Visible','off','units','normalized','position',[0.1,0.1,1,1]);
H1_hat5_tsne = tsne(H_hat5(:,1:K1));
subplot(1,2,1);
gscatter(H1_hat5_tsne(:,1),H1_hat5_tsne(:,2),cell_label);
title('tsne H1');


H12_hat5_tsne = tsne(H_hat5(:,1:K1+K2));
subplot(1,2,2);
gscatter(H12_hat5_tsne(:,1),H12_hat5_tsne(:,2),cell_label);
title('tsne H12');

saveas(gcf,sprintf('%s/result_H12.png',res_path));
% 
% figure('Visible','off','units','normalized','position',[0.1,0.1,1,1]);
% for k = 1:10
%     subplot(3,4,k);
%      gscatter(1:n,H_train(K1+k,:),cell_label);
% %     scatter(1:n,H_train(K1+k,:),100,cell_color,'.');
%     legend_s = findobj('type','legend');
%     delete(legend_s)
% 
%     title(sprintf('%d', k));
% end
% saveas(gcf,sprintf('%s/Check_H_110.png',res_path));
% 
% figure('Visible','off','units','normalized','position',[0.1,0.1,1,1]);
% for k = 1:10
%     subplot(3,4,k);
%     gscatter(1:n,H_train(K1+10+k,:),cell_label);
%         legend_s = findobj('type','legend');
%     delete(legend_s)
% 
% %     scatter(1:n,H_train(K1+10+k,:),100,cell_color,'.');
%     title(sprintf('%d', k+10));
% end
% saveas(gcf,sprintf('%s/Check_H_1120.png',res_path));


figure('Visible','off','units','normalized','position',[0.1,0.1,1,1]);
subplot(3,3,1)
plot(1:K1,W1_train_norm,'-oc');
legend('W1');
xlabel('Components');
ylabel('L2 norm');
title('W1');

subplot(3,3,2)
plot(1:K2,W2_train_norm,'-or');
legend('W2');
xlabel('Components');
ylabel('L2 norm');
title('W2');

subplot(3,3,3)
plot(1:K3,W3_train_norm,'-oy');
legend('W3');
xlabel('Components');
ylabel('L2 norm');
title('W3');


subplot(3,3,4)
plot(1:K1,H1_train_norm,'-oc');
legend('H1');
xlabel('Components');
ylabel('L2 norm');
title('H1');

subplot(3,3,5)
plot(1:K2,H2_train_norm,'-or');
legend('H2');
xlabel('Components');
ylabel('L2 norm');
title('H2');

subplot(3,3,6)
plot(1:K3,H3_train_norm,'-oy');
legend('H3');
xlabel('Components');
ylabel('L2 norm');
title('H3');

subplot(3,3,7)
plot(1:K2,alpha_H2_train,'-*b');
legend('alpha 2');
xlabel('Components');
ylabel('alpha value');
title('alpha 2');

subplot(3,3,8)
plot(1:K3,alpha_H3_train,'-*g');
legend('alpha 3');
xlabel('Components');
ylabel('alpha value');
title('alpha 3');

saveas(gcf,sprintf('%s/Check_W_A_n%dp%dsigma%0.7f.png',res_path,n,p,sigma_setting));
% fprintf('\nPLOT FINISHED\n');

end