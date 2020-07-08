clear;
clc;
K3 = 5;
multi_donor = false;
bulk_name = 'blood';
datanames = {'donor_BM0828'};
peak_selection_k = 3;

for data_name = datanames
for K2 = [5,10,15]
RA3_script(cell2mat(data_name), bulk_name, peak_selection_k, K2, K3, multi_donor);
end
end
