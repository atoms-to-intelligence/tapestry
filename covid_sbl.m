%% Initialization
clear;
clc;
matrices;
[m,n] = size (M5); A = M5;

% m = 40;
% n = 100;
% mr = ceil(0.9*m);
% mcv = m-mr;
for s = 1:10
nsignals = 100;

% A = rand(m,n); A(A < 0.5) = 0; A(A >= 0.5) = 1;

FN_count = 0;
FP_count = 0;
for kk = 1:nsignals
    x = zeros(n,1);
    indices = randperm(n); x(indices(1:s)) = 999*rand(s,1)+1; true_supp = indices(1:s);
    Ax = A*x;
    y = Ax;

%     sigval = 0.01*median(abs(Ax(Ax>0)));
%     y(y>0) = y(y>0) + randn(length(y(y>0)),1)*sigval;
    
    y = y + randn(m,1).*(0.01*y);
    sigval = 0.01*norm(y,2);
%     yr = y(1:mr);
%     ycv = y(mr+1:m);
%     Ar = A(1:mr,:);
%     Acv = A(mr+1:m,:);

    tau = 0.01*min(x(x>0)); % we would know the least virus concentration in an infected person. set tau to be 0.01 times that
                            % tau will be required to remove insignificant entries of x_est.

    %% Pre-Processing
    nz_index = find(y);
    z_index = find(~y);
    red_y = y(nz_index);
%     red_A = Ar(nz_index,:);

    [r,c] = find(A(z_index,:));
    ind_zero_x = unique(c);
    ind = [1:length(x)];
    ind_nonzero_x = setxor(ind,intersect(ind,ind_zero_x));
    red_x = x(ind_nonzero_x);
    red_A = A(nz_index,ind_nonzero_x);

    red_n = length(red_x);


    %% Sparse Bayesian Learning
    if red_n ==0
        x_est = zeros(red_n,1);
    else
        inv_Gamma = eye(red_n);

        mu_old = ones(red_n,1);
        mu = zeros(red_n,1);
        while norm(mu_old-mu,2)/norm(mu,2)>1e-4
            mu_old = mu;
            inv_Sigma = sigval^(-2)*(red_A'*red_A) + inv_Gamma;
            mu = sigval^(-2)*(inv_Sigma\(red_A'*red_y));
            Sigma = inv(inv_Sigma);
            gamma = mu.^2 + diag(Sigma);
            inv_Gamma = diag(1./gamma);
            rmse = norm(mu-red_x,2)/norm(red_x,2);
%             fprintf("\n Iter: %d, rmse= %f",iter,rmse)s
        end
        x_est = zeros(n,1);
        x_est(ind_nonzero_x) = mu;
        x_est(x_est < tau) = 0;
        
    end
        rmse(kk) = norm(x_est-x,2)/norm(x,2);
        FN = length(find(x_est == 0 & x > 0));
        FP = length(find(x_est > 0 & x == 0));
        if FN==0
            FN_count = FN_count+1;
        end
        if FP==0
            FP_count = FP_count+1;
        end
        false_neg(kk) = FN/s;
        false_pos(kk) = FP/(n-s);
        accuracy(kk) =  length(find((x_est > 0 & x > 0) |(x_est == 0 & x == 0) ))/n;
%         fprintf('\nRMSE: %f, FN: %f, FP: %f, Acc: %f',rmse(kk),false_neg(kk),false_pos(kk),accuracy(kk));
end

fprintf('\ns = %d, avg. RMSE = %f, avg. FNR = %f, avg. FPR = %f, avg. Acc = %f',s,mean(rmse),mean(false_neg),mean(false_pos),mean(accuracy));
fprintf('\n No FN in %d cases, No FP in %d cases ',FN_count,FP_count)
end