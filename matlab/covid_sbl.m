clear;
clc;
matrices;

%% Refer to SLide 31 of http://math.iisc.ernet.in/~nmi/Chandra%20Murthy.pdf for algorithm
%% Initialization
[m,n] = size (M2); A = M2;

% m = 40;
% n = 100;
% mr = ceil(0.9*m);
% mcv = m-mr;
for s = 1:10
nsignals = 1000;

% A = rand(m,n); A(A < 0.5) = 0; A(A >= 0.5) = 1;


fn = zeros(nsignals,1); fp = zeros(nsignals,1); tn = zeros(nsignals,1); tp = zeros(nsignals,1);
noFP_count = 0; noFN_count=0; noError_count=0;

for kk = 1:nsignals
    x = zeros(n,1);
    % Generate random signals,x with desired sparsity, s
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
    % As all A_ij and x_j are positive, for any y_i=0 implies that for all j s.t A_ij=1, x_j=0. 
    % This reduces problem dimension.
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
        % E-step
        %   mu is estimated mean of posterior distribution x|y, and so is the estimated red_x computed iteratively
        %   Sigma is variance of posterior distribution of x|y
        % M-step
        %   Gamma is the prior variance of x, inv_Gamma is saved as E-step only requires the inverse
        mu_old = ones(red_n,1);
        mu = zeros(red_n,1);
        while norm(mu_old-mu,2)/norm(mu,2)>1e-4
            mu_old = mu;
            inv_Sigma = sigval^(-2)*(red_A'*red_A) + inv_Gamma;
            mu = sigval^(-2)*(inv_Sigma\(red_A'*red_y));
            Sigma = inv(inv_Sigma);
            gamma = mu.^2 + diag(Sigma);
            inv_Gamma = diag(1./gamma);
            rmse = norm(mu-red_x,2)/norm(red_x,2); %RMS error of estimation compared to true signal
%             fprintf("\n Iter: %d, rmse= %f",iter,rmse)s
        end
        x_est = zeros(n,1);
        x_est(ind_nonzero_x) = mu;
        x_est(x_est < tau) = 0;
        
    end
        rmse(kk) = norm(x_est-x,2)/norm(x,2);
        fn(kk) = length(find(x_est == 0 & x > 0));
        fp(kk) = length(find(x_est > 0 & x == 0));
        tn(kk) = length(find(x_est == 0 & x == 0));
        tp(kk) = length(find(x_est > 0 & x > 0));
        if fn(kk) ==0
            noFN_count = noFN_count+1;
        end
        
        if fp(kk) ==0
            noFP_count = noFP_count+1;
        end
        if (fn(kk)+fp(kk)) ==0
            noError_count = noError_count+1;
        end
        
end
FN = sum(fn);
FP = sum(fp);
TN = sum(tn);
TP = sum(tp);
if isempty(fn(fn>1)) ==1 
    avg_fn = 0;
else
    avg_fn = mean(fn(fn>1));
end
if isempty(fp(fp>1)) ==1
    avg_fp = 0;
else
    avg_fp = mean(fp(fp>1));
end
precision = TP/(TP+FP);
recall = TP/(TP+FN);
% fprintf('\ns = %d, avg. RMSE = %f, avg. Precision = %f, avg. Recall = %f',s,mean(rmse),precision,recall);
fprintf(' %d,%f,%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d \n',s,precision,recall,avg_fp,avg_fn,FP,FN,TP,(TP+FN),(TP+FP),noFP_count,noFN_count, noError_count)
end
