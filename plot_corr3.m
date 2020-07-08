function new_score = calculatePCs(data,coeff,center)
% Author: Jason Buenrostro, Broad Institute
% calculates PC scores

% number of components
num_comp = size(coeff,2);

% center
if strcmp(center,'center')==1
    mean_data = mean(data);
    center_data = data - repmat(mean_data,[size(data,1) 1]);
else
    disp('Not centered!')
    center_data = data;
end

% calc score
clear new_score
for j = 1:num_comp
    % all data
    disp(j)
    for i = 1:size(data,2)
        new_score(i,j) = sum(coeff(:,j).*center_data(:,i));
    end
end

end








clc
clear
rng(23);

%%
data_name = 'donor_BM0828';
bulk_name = 'blood_RSHscMppLmppCmp';

res_path = sprintf('../resultForComparison/results/%s',datestr(now,30));
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

% train PCA on bulk data
[coeff,score,~,~,explained] = pca(bulk_mat,'Centered',true,'economy',true);
% [coeff,score,~,~,explained] = pca(cqnAll.newDat(:,distalIdx),'Centered',true,'economy',true);

%% score scATAC data
% calc score
data = Y;
% data = full(scDat.counts(:,distalIdx))';
sc_PCscores = calculatePCs(data,coeff,'center');

% pca
corrMat = 1-squareform(pdist(sc_PCscores,'correlation'));
[a,b,~,~,e] = pca(corrMat,'NumComponents',10);
sum(e(1)+e(2))

load('../../supplementary_code/data/donor_BM0828_cell.mat')
cell_color = donor_BM0828_cell;
%% plot!
% plot and change defaults to match publication plot
figure;scatter3(b(:,1),b(:,2),b(:,3),100,squeeze(cell_color),'.');
% figure;scatter3(b(:,1),b(:,2),b(:,3),100,squeeze(scDat.colors(:,:,:)),'.');
view(80,-15);set(gca,'Xdir','reverse','Ydir','reverse','Zdir','reverse')
xlabel('PC1');ylabel('PC2');zlabel('PC3')
set(gca, 'FontSize',15, 'LineWidth',1.5)
saveas(gcf,sprintf('%s/Corr_scatter3.png',res_path));

