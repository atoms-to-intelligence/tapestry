n = 64;
m = 24;
d = 3;
A = zeros(m,n);
x = 11;
for i = 1:m
    r = randperm(n,x);
    A(i,r) = ones(1,x);
end
%A = randi(2,m,n)-1;
num = 100000;

true_positives = 0;
predicted_positives = 0;
actual_positives = num*d;
tau = 0;
neg = 0;

for it = 1:num

    x_true = zeros(n,1);
    r = randperm(n,d);
    x_true(r) = randi(10000,d,1);
    
    y = A*x_true;
    
    residual = y;
    support = zeros(n);
    epsilon = 0.00001;
    
    while(norm(r)>epsilon)
        temp = norm(residual);
        [a,j] = max(abs( residual' * (A ./ vecnorm(A)) ));
        support(j) = 1;
        theta = pinv( A(:, support == 1) ) * y;
        residual = y - A(:, support == 1) * theta;
        if norm(residual) == temp
           break
       end
    end
    
    x_recon = zeros(n,1);
    x_recon(support == 1) = theta;
    
    predicted_positives = predicted_positives + sum(x_recon>tau);
    neg = neg + sum(x_recon<0);
    
    for i = 1:n
        if x_true(i) > 0 && x_recon(i) > tau
            true_positives = true_positives+1;
        end
    end
    
end

precision = true_positives / predicted_positives;
recall = true_positives / actual_positives;