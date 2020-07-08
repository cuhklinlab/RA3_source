function [W_o,H_o,Gamma_o,lgp_o,A_o,Beta_o,sigma_s_o] = RA3_func(donor_label,cell_label,Y,theta,tau,tau0,tau1,K1,K2,K3,Gamma,A,W,V,sigma_s,H,X,Beta,res_path,set_sparse,fix_A)
% profile on;

% tic;
fid = fopen(sprintf('%s/res.txt',res_path),'wt');
MAXIte = 1000;
[p,n] = size(Y);

K = K1 + K2 + K3;
err = 1e-6;
Qy = -inf(1,MAXIte);
Qw = -inf(1,MAXIte);
Qh = -inf(1,MAXIte);
Qgamma = -inf(1,MAXIte);
log_p = -inf(1,MAXIte);
diff_log_Q = -inf(1,MAXIte);
Q_H = -inf(1,MAXIte);
max_A = [1e12*ones(1,K2) 1e12*ones(1,K3)];
cont = zeros(1,n); % contribution of each sample for logp

phi_gamma_0 = tau0^-2;
phi_gamma_1 = tau1^-2;
for ite = 2:MAXIte
    
%% Update H
    E_h  = zeros(K,n);
    E_hh = zeros(K,K,n);
    W_s = 1/sigma_s * W' * W;
    D = (ones(K,n) - Gamma) * tau0^2 + Gamma * tau1^2;
    D(1:K1,:) = tau^2 * ones(K1,n);
    D(K1+K2+1:K,:) = tau^2 * ones(K3,n);
    for j = 1:n
        D_j = D(:,j);

        SIGMA_h_inv = W_s + diag(D_j.^(-1)); % K*K
        SIGMA_h = (SIGMA_h_inv)^(-1); % K*K
        MU_h = 1/sigma_s *(Y(:,j)-Beta*X(:,j))' * W * SIGMA_h; % 1*K
        E_h(:,j) = MU_h'; % K*1
        E_hh(:,:,j)  = E_h(:,j) * E_h(:,j)' + SIGMA_h; % K*K
    end

    
%     log_Q1 = Q_func(Y,Beta,X,W,E_h,E_hh,n,p,sigma_s,K1,K2,K3,K,Gamma,tau0,tau1,A,V,theta);
    
%% Update W
%     W = ((Y-Beta*X) * E_h' + sigma_s * V * A) / (sum(E_hh,3) + sigma_s * A);
    W(:,K1+1 : K) = ((Y-Beta*X-V(:,1:K1)*E_h(1:K1,:)) * E_h(K1+1 : K,:)' + sigma_s * V(:,K1+1 : K) * A) / (sum(E_hh(K1+1 : K, K1+1 : K,:),3) + sigma_s * A);

%% Update Beta
    Beta  = (Y-W*E_h) * X' / (X*X');    

%% Update A

    if(fix_A == false)
        A_new = diag(A);
        for k = 1 : K2+K3
            tmp = W(:,k+K1) - V(:,k+K1);
            A_new(k) = p / (tmp' * tmp);
            if A_new(k) > max_A(k)
                A_new(k) = max_A(k);
            end
            
        end
        A = diag(A_new);
    end
    
%% Update Gamma
    if(set_sparse == true)
        Gamma_new = ones(K,n);
        for k = (K1+1):(K1+K2)
            for j = 1:n
                f_gamma_0 = - 0.5 * (E_hh(k,k,j)) * phi_gamma_0 + 0.5 * log(phi_gamma_0) + log(1-theta);
                f_gamma_1 = - 0.5 * (E_hh(k,k,j)) * phi_gamma_1 + 0.5 * log(phi_gamma_1) + log(theta);
                Gamma_new(k,j) = 1 * (f_gamma_1 > f_gamma_0);
            end
        end
        Gamma = Gamma_new;
        % disp(sprintf('\t\tsum_Gamma: %s\tmean_Gamma: %s', num2str(sum(Gamma,'all')), num2str(mean(Gamma,'all')) ));
    end
    
%% Update sigma_square
    sigma_s = 0;
    W_tmp = W' * W;
    TMP = (Y - Beta*X); % try
    Mu_tmp = W' * TMP; % try
    for j = 1:n
        sigma_s = sigma_s + TMP(:,j)' * TMP(:,j) - 2*Mu_tmp(:,j)'*E_h(:,j) + trace(E_hh(:,:,j)*W_tmp);
    end
    sigma_s = sigma_s / (p*n); 
    % log_Q2 = Q_func(Y,Beta,X,W,E_h,E_hh,n,p,sigma_s,K1,K2,K3,K,Gamma,tau0,tau1,A,V,theta);
    
    
%% Calculate log_p

    Qw(ite) = - (K2+K3) * p/2 * log(2*pi) + p/2 * trace(log(A)) - 1/2*trace(A*(W(:,K1+1 : K)-V(:,K1+1 : K))'*(W(:,K1+1 : K)-V(:,K1+1 : K)));
    Qgamma(ite) = sum(sum(Gamma(K1+1:K1+K2,:)*log(theta) + (ones(K2,n) - Gamma(K1+1:K1+K2,:))*log(1-theta)));
    log_tmp = - n * p/2 * log(2*pi*sigma_s);
    W_s = 1/sigma_s * W_tmp;
    D = (ones(K,n) - Gamma) * tau0^2 + Gamma * tau1^2;
    D(1:K1,:) = tau * ones(K1,n);
    D(K1+K2+1:K,:) = tau * ones(K3,n);
        
    for j = 1:n
        D_j = D(:,j);
        Sigma_j_inv = W_s + diag(D_j.^(-1));
        Sigma_j = (Sigma_j_inv)^(-1);
        Mu_j = 1/sigma_s *  Sigma_j* Mu_tmp(:,j);       
        cont(j) = - 1/2*log(prod(D_j)) - 1/2*(1/sigma_s * TMP(:,j)' * TMP(:,j) ) + 1/2*log(det(Sigma_j)) + 1/2*Mu_j'* Sigma_j^(-1)*Mu_j; % 2nd
        log_tmp = log_tmp + cont(j);
    end
    
    Qh(ite) = log_tmp;
    log_p(ite) = log_tmp + Qw(ite) + Qgamma(ite);
    
    if(log_p(ite) > log_p(ite-1))
        H_o = E_h;
        W_o = W;
        Beta_o = Beta;
        Gamma_o = Gamma;
        A_o = A;
        sigma_s_o = sigma_s;
        lgp_o = log_p(ite);
    else
        fprintf(fid,'[Warning] Reduced log_p.');
        error('[Warning] Reduced log_p.');
    end
        
    fprintf(fid, '[%d] Time: %0.1fs\tH: %s\tW: %s\tGamma: %s\tlog_p: %s\n',ite,toc,num2str(Qh(ite)),num2str(Qw(ite)),num2str(Qgamma(ite)),num2str(log_p(ite))  );


    if(abs(log_p(ite) - log_p(ite-1)) < err * abs(log_p(ite-1)))
        break;
    end
end
fprintf(fid, '\nSetting\nN: %d\tP: %d\tK: %d\tTime: %0.1fs\tmax_log_p: %s\n',n,p,K,toc,num2str(max(log_p))  );

scale_W = sqrt(diag(W_o' * W_o))'; 
% size(scale_W)% 1 by K
W_o = W_o ./ scale_W;
H_o = H_o .* scale_W';

end