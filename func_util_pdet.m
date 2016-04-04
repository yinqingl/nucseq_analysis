%{
helper functions to be used with the probability of detection

Yinqing Li
yinqingl@mit.edu
2015

%import functions
fh = func_util_pdet;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

Copyright (c) 2015, Yinqing Li
%}
function func_h = func_util_pdet
func_h = struct();
func_h.pdist_w = @pdist_w;  %pairwise distance with weights between columns from input X, see pdist
func_h.pdist2_w = @pdist2_w;    %pairwise distance with weights between columns from inputs X and Y, see pdist2
func_h.corr_w = @corr_w;    %pairwise correlation with weights between columns from input X, see corr
func_h.costheta_w = @costheta_w; %pairwise cosine with weights between columns from input X, similar to corr
func_h.cov_w = @cov_w; %pairwise covariance with weights between columns from input X, see cov
func_h.cov_w_nz = @cov_w_nz;  %pairwise inner product with weights between columns from input X
func_h.zeromean_w = @zeromean_w; %substract columns of X by their mean calculated with weights
func_h.pca_w = @pca_w; %PCA with weights
func_h.ttest2_w = @ttest2_w; %t test with weights, see ttest2
func_h.quantile_w = @quantile_w; %quantiles with weights, see quantile
func_h.mean_w = @mean_w; %mean with weights, see mean
func_h.std_w = @std_w; %standard deviation with weights, see std
func_h.std_robust_w = @std_robust_w; %standard deviation with weights using Rousseeuw and Croux method, see std
func_h.median_w = @median_w; %median with weights, see median
func_h.pmean_w = @pmean_w; %pairwise average of elements of a vector
func_h.prepEcEr_isa2_w = @prepEcEr_isa2_w;  %prepare Ec and Er for isa2
% func_h.spatialcrosscorr_w = @spatialcrosscorr_w;
func_h.spatialcorr_w = @spatialcorr_w;  %crosscorrelation with weights and spatial distance on the tsne embedding
func_h.lratiotest2_binorm_w = @lratiotest2_binorm_w; %two-sample likelihood ratio test using combined binomial and normal model
func_h.zttest2_w = @zttest2_w; %two-sample binomial and normal test

%make sure these functions are in the directory
%{
estimate_detection_probability_kde; %estimate pdet given expression data
%}

end

function example_code

option = struct();
option.is_show_iteration = true;
option.max_iter = 300;
option.beta_lb = [-1e-4,1e-4];
option.alg = 'geometric';

[dat_p_det,dat_prob_est] = estimate_detection_probability_kde(dat,option);
W = ones(size(dat));
W(dat == 0) = dat_p_det;
lambda_ik = dat_prob_est.lambda; %detection probability function
beta_i = dat_prob_est.beta_i; %parameters for detection probability function

%plot cell detection probability
figure
for i = [1:10],
x = exp(linspace(0,6,100));
plot(log(x+1),lambda_ik(x',beta_i(:,i))); hold on;
end

end

%{
pairwise distance with weights between columns from input X, see pdist

inputs:
if a is given, then calculate pdist on y
defined by linear transformation y = x*a;
rows of x correspond to observations (genes)
%}
function d = pdist_w(X,W,DISTANCE,a)
X = double(X); W = double(W);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end
if strcmp(DISTANCE,'euclidean'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);
        
%         mode = 3;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        [n,p] = size(X);
        [~,m] = size(a);
        
        tic
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:n-1,
            dsq = zeros(n-i,m);
            w_sum = zeros(n-i,1);
            x_i = X(i,:);
            y_i = X((i+1):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i+1):n,:);

            for q = 1:p     %jit
                w_j = w_x_i(q)*w_y_i(:,q);
                d_xy = (x_i(q) - y_i(:,q)).*w_j;
                for kk = 1:m,
                     dsq(:,kk) = dsq(:,kk) + d_xy.*a(q,kk);
                end
                w_sum = w_sum + w_j;
            end
            w_sum = w_sum/p;
            for kk = 1:m,
                dsq(:,kk) = dsq(:,kk)./w_sum;
                dsq(:,kk) = dsq(:,kk).*dsq(:,kk);
            end    
            dsq = sum(dsq,2);
            d(k:(k+n-i-1)) = dsq;
            k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end            
        end
        d = sqrt(d);
%         toc

    else
        [n,p] = size(X);

        tic
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:n-1,
            dsq = zeros(n-i,1);
            w_sum = zeros(n-i,1);
            x_i = X(i,:);
            y_i = X((i+1):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i+1):n,:);    
            for q = 1:p     %jit
                w_j = w_x_i(q)*w_y_i(:,q);
                d_xy = x_i(q) - y_i(:,q);
                dsq = dsq + d_xy.*d_xy.*w_j;
                w_sum = w_sum + w_j;
            end
            d(k:(k+n-i-1)) = dsq./w_sum;
            k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
        d = sqrt(d*p);
%         toc      
    end
elseif strcmp(DISTANCE,'correlation'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);
%         mode = 4;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        d = corr_w(X', W', a');
    else
%         mode = 1;
%         d = pdet_util_c(X,W,mode);
%         d = squareform(d);

        d = corr_w(X', W');
    end
    d = squareform(1-d);
elseif strcmp(DISTANCE,'cosine'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);         
%         mode = 4;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        d = costheta_w(X', W', a');
    else
%         mode = 1;
%         d = pdet_util_c(X,W,mode);
%         d = squareform(d);

        d = costheta_w(X', W');
    end
    d = squareform(1-d);
else
    disp('missing DISTANCE');
end

end

%{
pairwise distance with weights between columns from inputs X and Y, see pdist2

inputs:
if a is given, then calculate pdist on zx and zy
defined by linear transformation zx = x*a, zy = y*a;
rows of x,y correspond to observations (genes)
%}
function d = pdist2_w(X,Y,W_X,W_Y,DISTANCE,a)
X = double(X); Y = double(Y);
W_X = double(W_X); W_Y = double(W_Y);
% is_a = 1;
% if ~exist('a', 'var')
%     is_a = 0;
% end
% x = [X;Y];
% w = [W_X;W_Y];
% if is_a,
%     a = double(a);
%     d = squareform(pdist_w(x,w,DISTANCE,a));
% else,
%     d = squareform(pdist_w(x,w,DISTANCE));
% end
% d = d(1:size(X,1),size(X,1)+1:end);

is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end
if strcmp(DISTANCE,'euclidean'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);
        
%         mode = 3;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        [n,p] = size(X);
        [ny,py] = size(Y);
        [~,m] = size(a);
        
        tic
        d = zeros(n,ny);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
%         for i = 1:n-1,
            dsq = zeros(ny,m);
            w_sum = zeros(ny,1);
            x_i = X(i,:);
            y_i = Y;
            w_x_i = W_X(i,:);
            w_y_i = W_Y;

            for q = 1:p     %jit
                w_j = w_x_i(q)*w_y_i(:,q);
                d_xy = (x_i(q) - y_i(:,q)).*w_j;
                for kk = 1:m,
                     dsq(:,kk) = dsq(:,kk) + d_xy.*a(q,kk);
                end
                w_sum = w_sum + w_j;
            end
            w_sum = w_sum/p;
            for kk = 1:m,
                dsq(:,kk) = dsq(:,kk)./w_sum;
                dsq(:,kk) = dsq(:,kk).*dsq(:,kk);
            end    
            dsq = sum(dsq,2);
            d(i,:) = dsq;
%             d(k:(k+n-i-1)) = dsq;
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end            
        end
        d = sqrt(d);
%         toc

    else
%         mode = 0;
%         d = pdet_util_c(X,W,mode);
%         d = squareform(d);
        
        [n,p] = size(X);
        [ny,py] = size(Y);
        
        tic
        d = zeros(n,ny);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
%         for i = 1:n-1,
            dsq = zeros(ny,1);
            w_sum = zeros(ny,1);
            x_i = X(i,:);
            y_i = Y;
            w_x_i = W_X(i,:);
            w_y_i = W_Y;    
            for q = 1:p     %jit
                w_j = w_x_i(q)*w_y_i(:,q);
                d_xy = x_i(q) - y_i(:,q);
                dsq = dsq + d_xy.*d_xy.*w_j;
                w_sum = w_sum + w_j;
            end
            d(i,:) = dsq./w_sum;
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
        d = sqrt(d*p);
%         toc      
    end
elseif strcmp(DISTANCE,'correlation'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);
%         mode = 4;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        d = corr2_w(X', W_X', Y', W_Y', a');
        
    else
%         mode = 1;
%         d = pdet_util_c(X,W,mode);
%         d = squareform(d);

        d = corr2_w(X', W_X', Y', W_Y');
    end
    d = 1-d;
%     d = squareform(1-d);
elseif strcmp(DISTANCE,'cosine'),
    if (is_a),
        %assuming this is from PCA
        %X should be gene centered
%         X = zeromean_w(X, W);         
%         mode = 4;
%         d = pdet_util_c(X,W,mode,a);
%         d = squareform(d);
        
        d = costheta2_w(X', W_X', Y', W_Y', a');
    else
%         mode = 1;
%         d = pdet_util_c(X,W,mode);
%         d = squareform(d);

        d = costheta2_w(X', W_X', Y', W_Y');
    end
    d = 1-d;
%     d = squareform(1-d);
else
    disp('missing DISTANCE');
end

end

%{
pairwise correlation with weights between columns from input X, see corr

inputs:
if a is given, then calculate pdist on y
defined by linear transformation

y = X'*a'; X should be gene centered, correlation calculated on columns
if y = X*a', then correlation calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = corr_w(X, W, a)
X = double(X); W = double(W);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W = W';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [~,m] = size(a);
        
        X = zeromean_w(X, W);
        
        tic
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:n-1,
            xw = zeros(n-i,m);
            yw = zeros(n-i,m);
            xy = zeros(n-i,1);
            xx = zeros(n-i,1);
            yy = zeros(n-i,1);

            u_x = zeros(n-i,1);
            u_y = zeros(n-i,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = X((i+1):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i+1):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

            for kk = 1:m,
                u_x = u_x + xw(:,kk);
                u_y = u_y + yw(:,kk);
            end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
            u_x = u_x/m;
            u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
                xx = xx + xw(:,kk).*xw(:,kk);
                yy = yy + yw(:,kk).*yw(:,kk);                
            end
            xy = xy - u_x.*u_y*m;
            xx = xx - u_x.*u_x*m;
            yy = yy - u_y.*u_y*m;
            
            d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
            k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
%         toc
        rho = squareform(d);
        rho(1:n+1:n*n) = 1.0;
    else
        %approximate center
        %exact center would use pairwise w_j
        X = zeromean_w(X, W);
        X = X.*W;
        y = X*a';
        rho = corr(y);
    end
else
    X = X';
    W = W';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    
    tic
    d = zeros(1,n*(n-1)./2);
    k = 1;
    for i = 1:n-1,
        xy = zeros(n-i,1);
        xx = zeros(n-i,1);
        yy = zeros(n-i,1);

        u_x = zeros(n-i,1);
        u_y = zeros(n-i,1);

        w_sum = zeros(n-i,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = X((i+1):n,:);
        w_x_i = W(i,:);
        w_y_i = W((i+1):n,:);

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
            w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
            u_x = u_x + xw;
            u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
            xx = xx + x_i(q).*xw;
            yy = yy + y_i(:,q).*yw;            
        end
        u_x = u_x./w_sum;
        u_y = u_y./w_sum;
        xy = xy - u_x.*u_y.*w_sum;
        xx = xx - u_x.*u_x.*w_sum;
        yy = yy - u_y.*u_y.*w_sum;

        d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
        k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end            
    end

    rho = squareform(d);
    rho(1:n+1:n*n) = 1.0;
    
end

end

%{
pairwise correlation with weights between columns from inputs X and Y, see corr2

inputs:
if a is given, then calculate pdist on z
defined by linear transformation

%zx = X'*a'; X should be gene centered, correlation calculated on columns
%if zx = X*a', then correlation calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = corr2_w(X, W_X, Y, W_Y, a)
X = double(X); Y = double(Y);
W_X = double(W_X); W_Y = double(W_Y);

is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W_X = W_X';
        Y = Y';
        W_Y = W_Y';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [ny,py] = size(Y);
        [~,m] = size(a);
        
%         X = zeromean_w(X, W_X);
        
        tic
        d = zeros(n,ny);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
%         for i = 1:n-1,
            xw = zeros(ny,m);
            yw = zeros(ny,m);
            xy = zeros(ny,1);
            xx = zeros(ny,1);
            yy = zeros(ny,1);

            u_x = zeros(ny,1);
            u_y = zeros(ny,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = Y;
%             y_i = X((i+1):n,:);
            w_x_i = W_X(i,:);
            w_y_i = W_Y;
%             w_y_i = W((i+1):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

            for kk = 1:m,
                u_x = u_x + xw(:,kk);
                u_y = u_y + yw(:,kk);
            end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
            u_x = u_x/m;
            u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
                xx = xx + xw(:,kk).*xw(:,kk);
                yy = yy + yw(:,kk).*yw(:,kk);                
            end
            xy = xy - u_x.*u_y*m;
            xx = xx - u_x.*u_x*m;
            yy = yy - u_y.*u_y*m;
            
            d(i,:) = xy./(sqrt(xx.*yy)+realmin);
%             d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
%         toc
        rho = d;
%         rho = squareform(d);
%         rho(1:n+1:n*n) = 1.0;
    else
        %approximate center
        %exact center would use pairwise w_j
        X = zeromean_w(X, W_X);
        X = X.*W_X;
        Y = zeromean_w(Y, W_Y);
        Y = Y.*W_Y;
        
        x = X*a';
        y = Y*a';
%         rho = corr(y);
        rho = corr(x,y);
    end
else
    X = X';
    W_X = W_X';
    Y = Y';
    W_Y = W_Y';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    [ny,py] = size(Y);
    
    tic
    d = zeros(n,ny);
%     d = zeros(1,n*(n-1)./2);
%     k = 1;
    for i = 1:n,
%     for i = 1:n-1,
        xy = zeros(ny,1);
        xx = zeros(ny,1);
        yy = zeros(ny,1);

        u_x = zeros(ny,1);
        u_y = zeros(ny,1);

        w_sum = zeros(ny,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = Y;
%         y_i = X((i+1):n,:);
        w_x_i = W_X(i,:);
        w_y_i = W_Y;
%         w_y_i = W((i+1):n,:);

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
            w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
            u_x = u_x + xw;
            u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
            xx = xx + x_i(q).*xw;
            yy = yy + y_i(:,q).*yw;            
        end
        u_x = u_x./w_sum;
        u_y = u_y./w_sum;
        xy = xy - u_x.*u_y.*w_sum;
        xx = xx - u_x.*u_x.*w_sum;
        yy = yy - u_y.*u_y.*w_sum;
        
        d(i,:) = xy./(sqrt(xx.*yy)+realmin);
%         d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%         k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end            
    end
    
    rho = d;
%     rho = squareform(d);
%     rho(1:n+1:n*n) = 1.0;
    
end

end

%{
pairwise cosine with weights between columns from input X

inputs:
if a is given, then calculate pdist on y
defined by linear transformation

%y = X'*a'; X should be gene centered, cosine calculated on columns
%if y = X*a', then cosine calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = costheta_w(X, W, a)
X = double(X); W = double(W);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W = W';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [~,m] = size(a);
        
        tic
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:n-1,
            xw = zeros(n-i,m);
            yw = zeros(n-i,m);
            xy = zeros(n-i,1);
            xx = zeros(n-i,1);
            yy = zeros(n-i,1);

%             u_x = zeros(n-i,1);
%             u_y = zeros(n-i,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = X((i+1):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i+1):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

%             for kk = 1:m,
%                 u_x = u_x + xw(:,kk);
%                 u_y = u_y + yw(:,kk);
%             end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
%             u_x = u_x/m;
%             u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
                xx = xx + xw(:,kk).*xw(:,kk);
                yy = yy + yw(:,kk).*yw(:,kk);                
            end
%             xy = xy - u_x.*u_y*m;
%             xx = xx - u_x.*u_x*m;
%             yy = yy - u_y.*u_y*m;
            
            d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
            k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
%         toc
        rho = squareform(d);
        rho(1:n+1:n*n) = 1.0;
    else
        %weight has no effect if not centered
        X = X.*W;
        y = X*a';
        
        [n,p] = size(y);
        
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:p-1,
            y_sum = 0;
            y2_sum = 0;
            yy2_sum = 0;            
            for j = 1:n,
                y_sum = y_sum + y(j,i).*y(j,i+1:p);
                y2_sum = y2_sum + y(j,i)*y(j,i);
                yy2_sum = yy2_sum + y(j,i+1:p).*y(j,i+1:p);
            end
            d(k:(k+n-i-1)) = y_sum./(sqrt(y2_sum.*yy2_sum)+realmin);
            k = k + (n-i);
        end
        rho = squareform(d);
        rho(1:n+1:n*n) = 1.0;
    end        
else
    X = X';
    W = W';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    
    tic
    d = zeros(1,n*(n-1)./2);
    k = 1;
    for i = 1:n-1,
        xy = zeros(n-i,1);
        xx = zeros(n-i,1);
        yy = zeros(n-i,1);

%         u_x = zeros(n-i,1);
%         u_y = zeros(n-i,1);

%         w_sum = zeros(n-i,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = X((i+1):n,:);
        w_x_i = W(i,:);
        w_y_i = W((i+1):n,:);

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
%             w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
%             u_x = u_x + xw;
%             u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
            xx = xx + x_i(q).*xw;
            yy = yy + y_i(:,q).*yw;            
        end
%         u_x = u_x./w_sum;
%         u_y = u_y./w_sum;
%         xy = xy - u_x.*u_y.*w_sum;
%         xx = xx - u_x.*u_x.*w_sum;
%         yy = yy - u_y.*u_y.*w_sum;

        d(k:(k+n-i-1)) = xy./sqrt(xx.*yy);
        k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end            
    end

    rho = squareform(d);
    rho(1:n+1:n*n) = 1.0;
    
end

end

%{
pairwise cosine with weights between columns from inputs X and Y

inputs:
if a is given, then calculate pdist on y
defined by linear transformation

%y = X'*a'; X should be gene centered, cosine calculated on columns
%if y = X*a', then cosine calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = costheta2_w(X, W_X, Y, W_Y, a)
X = double(X); W_X = double(W_X);
Y = double(Y); W_Y = double(W_Y);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W_X = W_X';
        Y = Y';
        W_Y = W_Y';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [ny,py] = size(Y);
        [~,m] = size(a);
        
        tic
        d = zeros(n,ny);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
            xw = zeros(ny,m);
            yw = zeros(ny,m);
            xy = zeros(ny,1);
            xx = zeros(ny,1);
            yy = zeros(ny,1);

%             u_x = zeros(n-i,1);
%             u_y = zeros(n-i,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = Y;
            w_x_i = W_X(i,:);
            w_y_i = W_Y;

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

%             for kk = 1:m,
%                 u_x = u_x + xw(:,kk);
%                 u_y = u_y + yw(:,kk);
%             end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
%             u_x = u_x/m;
%             u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
                xx = xx + xw(:,kk).*xw(:,kk);
                yy = yy + yw(:,kk).*yw(:,kk);                
            end
%             xy = xy - u_x.*u_y*m;
%             xx = xx - u_x.*u_x*m;
%             yy = yy - u_y.*u_y*m;
            d(i,:) = xy./(sqrt(xx.*yy)+realmin);
%             d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
%         toc
        rho = d;
%         rho = squareform(d);
%         rho(1:n+1:n*n) = 1.0;
    else
        %weight has no effect if not centered
        disp('check input.');
        %{
        X = X.*W_X;
        Y = Y.*W_Y;
        x = X*a';
        y = Y*a';
        
        [n,p] = size(x);
        [ny,py] = size(y);
        
        d = zeros(1,n*(n-1)./2);
        k = 1;
        for i = 1:p-1,
            y_sum = 0;
            y2_sum = 0;
            yy2_sum = 0;            
            for j = 1:n,
                y_sum = y_sum + y(j,i).*y(j,i+1:p);
                y2_sum = y2_sum + y(j,i)*y(j,i);
                yy2_sum = yy2_sum + y(j,i+1:p).*y(j,i+1:p);
            end
            d(k:(k+n-i-1)) = y_sum./(sqrt(y2_sum.*yy2_sum)+realmin);
            k = k + (n-i);
        end
        rho = squareform(d);
        rho(1:n+1:n*n) = 1.0;
        %}
        rho = 0;
    end        
else
    X = X';
    W_X = W_X';
    Y = Y';
    W_Y = W_Y';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    [ny,py] = size(Y);
    
    tic
    d = zeros(n,ny);
%     d = zeros(1,n*(n-1)./2);
%     k = 1;
    for i = 1:n,
        xy = zeros(ny,1);
        xx = zeros(ny,1);
        yy = zeros(ny,1);

%         u_x = zeros(n-i,1);
%         u_y = zeros(n-i,1);

%         w_sum = zeros(n-i,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = Y;
        w_x_i = W_X(i,:);
        w_y_i = W_Y;

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
%             w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
%             u_x = u_x + xw;
%             u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
            xx = xx + x_i(q).*xw;
            yy = yy + y_i(:,q).*yw;            
        end
%         u_x = u_x./w_sum;
%         u_y = u_y./w_sum;
%         xy = xy - u_x.*u_y.*w_sum;
%         xx = xx - u_x.*u_x.*w_sum;
%         yy = yy - u_y.*u_y.*w_sum;
        d(i,:) = xy./sqrt(xx.*yy);
%         d(k:(k+n-i-1)) = xy./sqrt(xx.*yy);
%         k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end            
    end
    
    rho = d;
%     rho = squareform(d);
%     rho(1:n+1:n*n) = 1.0;
    
end

end

% %covariance calculated on columns
% %y = X'*a';
% %if y = X*a', then covariance calculated on rows
% function c = cov_w(X, W, a)
% X = double(X); W = double(W);
% is_a = 1;
% if ~exist('a', 'var')
%     is_a = 0;
% end
% 
% if (is_a),
%     a = double(a);
%     if size(X,1) == size(a, 2),
%         X = X';
%         W = W';
%         a = a';        
%         mode = 5;
%         c = pdet_util_c(X,W,mode,a);        
%     else
%         %approximate center
%         %exact center would use pairwise w_j
%         X = zeromean_w(X, W);
%         X = X.*W;
%         y = X*a';
%         c = cov(y);
%         z = sum(sum(W,2).^2);
%         n = sum(sum(ones(size(W)),2).^2);
%         c = c/z*n;
%     end
% else
%     X = X';
%     W = W';
%     mode = 2;
%     c = pdet_util_c(X,W,mode);
% end
% 
% end

%{
pairwise covariance with weights between columns from input X, see cov

inputs:
if a is given, then calculate pdist on y
defined by linear transformation

%y = X'*a'; X should be gene centered, covariance calculated on columns
%if y = X*a', then covariance calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = cov_w(X, W, a)
X = double(X); W = double(W);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W = W';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [~,m] = size(a);
        
        X = zeromean_w(X, W);
        
        tic
        d = zeros(n,n);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
            xw = zeros(n-i+1,m);
            yw = zeros(n-i+1,m);
            xy = zeros(n-i+1,1);
%             xx = zeros(n-i,1);
%             yy = zeros(n-i,1);

            u_x = zeros(n-i+1,1);
            u_y = zeros(n-i+1,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = X((i):n,:);
%             y_i = X((i+1):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i):n,:);
%             w_y_i = W((i+1):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

            for kk = 1:m,
                u_x = u_x + xw(:,kk);
                u_y = u_y + yw(:,kk);
            end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
            u_x = u_x/m;
            u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
%                 xx = xx + xw(:,kk).*xw(:,kk);
%                 yy = yy + yw(:,kk).*yw(:,kk);                
            end
            xy = xy - u_x.*u_y*m;
%             xx = xx - u_x.*u_x*m;
%             yy = yy - u_y.*u_y*m;
            
            d(i,i:end) = xy;
%             d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
        toc
        rho = d;
        rho = rho + rho';
        rho(1:n+1:end) = rho(1:n+1:end)/2;        
%         rho = squareform(d);
%         rho(1:n+1:n*n) = 1.0;
    else
        %approximate center
        %exact center would use pairwise w_j
        X = zeromean_w(X, W);
        X = X.*W;
        y = X*a';
        rho = cov(y);
        z = sum(sum(W,2).^2);
        n = sum(sum(ones(size(W)),2).^2);
        rho = rho/z*n;
%         rho = corr(y);
    end
else
    X = X';
    W = W';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    
    tic
    d = zeros(n,n);
%     d = zeros(1,n*(n-1)./2);
%     k = 1;
    for i = 1:n,
        xy = zeros(n-i+1,1);
%         xx = zeros(n-i,1);
%         yy = zeros(n-i,1);

        u_x = zeros(n-i+1,1);
        u_y = zeros(n-i+1,1);

        w_sum = zeros(n-i+1,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = X((i):n,:);
        w_x_i = W(i,:);
        w_y_i = W((i):n,:);

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
            w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
            u_x = u_x + xw;
            u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
%             xx = xx + x_i(q).*xw;
%             yy = yy + y_i(:,q).*yw;            
        end
        u_x = u_x./w_sum;
        u_y = u_y./w_sum;
        xy = xy - u_x.*u_y.*w_sum;
%         xx = xx - u_x.*u_x.*w_sum;
%         yy = yy - u_y.*u_y.*w_sum;
        
        d(i,i:end) = xy;
%         d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%         k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end            
    end
    rho = d;
    rho = rho + rho';
    rho(1:n+1:end) = rho(1:n+1:end)/2;
%     rho = squareform(d);
%     rho(1:n+1:n*n) = 1.0;
    
end


end

%{
pairwise inner product with weights between columns from input X

inputs:
if a is given, then calculate pdist on y
defined by linear transformation

%y = X'*a'; X should be gene centered, inner product calculated on columns
%if y = X*a', then inner product calculated on rows

rows of x,y correspond to observations (genes)
%}
function rho = cov_w_nz(X, W, a)
X = double(X); W = double(W);
is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
    if size(X,1) == size(a, 2),
        X = X';
        W = W';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center
        [n,p] = size(X);
        [~,m] = size(a);
        
%         X = zeromean_w(X, W);
        
        tic
        d = zeros(n,n);
%         d = zeros(1,n*(n-1)./2);
%         k = 1;
        for i = 1:n,
            xw = zeros(n-i+1,m);
            yw = zeros(n-i+1,m);
            xy = zeros(n-i+1,1);
%             xx = zeros(n-i,1);
%             yy = zeros(n-i,1);

%             u_x = zeros(n-i,1);
%             u_y = zeros(n-i,1);

%             w_sum = zeros(n-i,1);

            %speed up indexing
            x_i = X(i,:);
            y_i = X((i):n,:);
            w_x_i = W(i,:);
            w_y_i = W((i):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q)*w_y_i(:,q);
%                 w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
                y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
                    yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

%             for kk = 1:m,
%                 u_x = u_x + xw(:,kk);
%                 u_y = u_y + yw(:,kk);
%             end
%             w_sum = w_sum/p;
%             u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
%             u_x = u_x/m;
%             u_y = u_y/m;
            for kk = 1:m,
%                 x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
%                 xx = xx + x.*x;
%                 yy = yy + y.*y;
                xy = xy + xw(:,kk).*yw(:,kk);
%                 xx = xx + xw(:,kk).*xw(:,kk);
%                 yy = yy + yw(:,kk).*yw(:,kk);                
            end
%             xy = xy - u_x.*u_y*m;
%             xx = xx - u_x.*u_x*m;
%             yy = yy - u_y.*u_y*m;
            
            d(i,i:end) = xy;
%             d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%             k = k + (n-i);
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
%         toc
        rho = rho + rho';
        rho(1:n+1:end) = rho(1:n+1:end)/2;
%         rho = squareform(d);
%         rho(1:n+1:n*n) = 1.0;
    else
        %approximate center
        %exact center would use pairwise w_j
%         X = zeromean_w(X, W);
        X = X.*W;
        y = X*a';
        rho = y'*y/(size(y,1)); %biased estimator
        z = sum(sum(W,2).^2);
        n = sum(sum(ones(size(W)),2).^2);
        rho = rho/z*n;        
%         rho = cov(y);
%         rho = corr(y);
    end
else
    X = X';
    W = W';
%     mode = 1;
%     rho = pdet_util_c(X,W,mode);
    
    [n,p] = size(X);
    
    tic
    d = zeros(n,n);
%     d = zeros(1,n*(n-1)./2);
%     k = 1;
    for i = 1:n,
        xy = zeros(n-i+1,1);
%         xx = zeros(n-i,1);
%         yy = zeros(n-i,1);

%         u_x = zeros(n-i,1);
%         u_y = zeros(n-i,1);

        w_sum = zeros(n-i+1,1);

        %speed up indexing
        x_i = X(i,:);
        y_i = X((i):n,:);
        w_x_i = W(i,:);
        w_y_i = W((i):n,:);

        for q = 1:p,
            w_j = w_x_i(q)*w_y_i(:,q);
            w_sum = w_sum + w_j;
            xw = x_i(q).*w_j;
            yw = y_i(:,q).*w_j;
%             u_x = u_x + xw;
%             u_y = u_y + yw;
            xy = xy + x_i(q).*yw;
%             xx = xx + x_i(q).*xw;
%             yy = yy + y_i(:,q).*yw;            
        end
%         u_x = u_x./w_sum;
%         u_y = u_y./w_sum;
%         xy = xy - u_x.*u_y.*w_sum;
%         xx = xx - u_x.*u_x.*w_sum;
%         yy = yy - u_y.*u_y.*w_sum;
        
        d(i,i:end) = xy;
%         d(k:(k+n-i-1)) = xy./(sqrt(xx.*yy)+realmin);
%         k = k + (n-i);
        if rem(i,100) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end
    end
    rho = rho + rho';
    rho(1:n+1:end) = rho(1:n+1:end)/2;
%     rho = squareform(d);
%     rho(1:n+1:n*n) = 1.0;
    
end

end

%{
substract columns of X by the mean of the columns

rows of x corresponds to observations
columns of x corresponds to variables
%}
function X_c = zeromean_w(X, W)
m_X = sum(X.*W,1)./sum(W,1);
X_c = X - repmat(m_X,size(X,1),1);
end

%{
PCA with weights

inputs:
X: data matrix, row: samples, column: variables
W: weight matrix
Centered: true or false, use centered covariance
NumComponents: number of PCs to calculate

outputs:
coeff: data matrix, row: samples, column: principal components
score: linear combination of variables to principal components
latent: variance of principal components
%}
function [coeff, score, latent] = pca_w(X, W, varargin)
X = double(X); W = double(W);

paramNames = {'Centered','NumComponents'};
defaults   = {true, 0};
[vCentered, vNumComponents] = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%rows of x corresponds to observations
%{
x=rand(6,3);
x = x';

DOF = size(x,1)-1;
x_cc = x - repmat(mean(x,1),size(x,1),1);
[U,S,V] = svd(x_cc'*x_cc/DOF,0);   %svd on covariance matrix
pcvars = S;
pc = V;
zscores = x_cc*pc;

[pc, zscores, pcvars] = pca(x);
[coeff, score, latent, tsquare] = princomp(x); %same as pca
%}

%rows of X correspond to observations
%first find pc without using w


[pc, ~, pc_v] = pca(X, 'Economy', true, 'Centered', vCentered);

a = pc;
if vCentered,
    c = cov_w(X, W, a');
else
    c = cov_w_nz(X, W, a');
end
[~,S,V] = svd(c,0);
V = V(1:size(pc,2),:);
S = S(1:size(pc,2),1:size(pc,2));

latent = diag(S);
coeff = pc*V;
score = X*coeff;

%robust pca
if vNumComponents > 1,
    a = coeff;
    Xa = X*a;
    Xaw = bsxfun(@rdivide,Xa,size(W,2)./sum(W,2));
    addpath('RegL1-ALM')
    
    alg = 1;
    
    switch alg
        case 1
            %Zheng 2012
            [L U_est V_est L1_error] = RobustApproximation_M_UV_TraceNormReg(...
                Xaw',ones(size(Xaw')),vNumComponents,1e-3,1.5,1,0);

        case 2
            %Lin 2009
            [L, S] = inexact_alm_rpca(Xaw');
            
        case 3
            %Shu 2014
            addpath('lrslibrary-master\algorithms\lrr\ROSL')
            tol = 1e-5;
            maxIter = 100;
            lambda = 1/sqrt(size(Xaw,1)); %2e-3;
            [D,alpha,~,~] = inexact_alm_rosl(Xaw',vNumComponents,lambda,tol,maxIter);
            L = D*alpha;

        case 4
            %Xu 2012
            addpath('lrslibrary-master\algorithms\rpca\OP-RPCA')
            lambda = 0.35; % [0,1] 0.5(+lowrank) 0.4 0.3 0.2(+sparse)
            [L,S] = mr_pca_part(Xaw', ones(size(Xaw')), lambda);
    end

    [pc, ~, pcvars] = pca(L');
    coeff = coeff*pc;
    score = Xaw*pc;
    latent = pcvars(1:vNumComponents);
end
end

%{
two-sample t-test with weights
test each row across columns
rows of x corresponds to observations
%}
function [H,P] = ttest2_w(X,Y,W_X,W_Y)
X = double(X); Y = double(Y);
W_X = double(W_X); W_Y = double(W_Y);
nx = sum(W_X,1);
ny = sum(W_Y,1);
ux = sum(X.*W_X,1)./nx;
uy = sum(Y.*W_Y,1)./ny;
x = bsxfun(@minus,X,ux);
y = bsxfun(@minus,Y,uy);
sx = sum(x.^2.*W_X)./(nx-1);
% sx = sum(x.^2.*W_X)./(nx);
sy = sum(y.^2.*W_Y)./(ny-1);
% sy = sum(y.^2.*W_Y)./(ny);
t = (ux-uy)./sqrt(sx./nx+sy./ny);
df = (sx./nx + sy./ny).^2./((sx./nx).^2./(nx-1)+(sy./ny).^2./(ny-1));
% df = (sx./nx + sy./ny).^2./((sx./nx).^2./(nx)+(sy./ny).^2./(ny));
P = 2*tcdf(-abs(t),df);
% P = min(1-P,P);
H = zeros(size(P));
H(P<0.05) = 1;
end

%{
quantile with weights, see quantile
for each row, calculate quantile of columns
%}
function m = quantile_w(X,W,q)
X = double(X);
W = double(W);
[~,p] = size(X);
[Y,I] = sort(X,1);
W = bsxfun(@rdivide, W, sum(W,1));
for i = 1:p,
    W(:,i) = W(I(:,i),i);
end
W_sum = cumsum(W,1);
m = zeros(1,p);
for i = 1:p,
    m(i) = Y(find(W_sum(:,i)>q,1),i);
end
end

%{
median with weights, see median
for each row, calculate median of columns
%}
function m = median_w(X,W)
m = quantile_w(X,W,0.5);
end

%{
mean with weights, see mean
for each row, calculate mean of columns
%}
function m = mean_w(X,W)
X = double(X);
W = double(W);
m = sum(X.*W,1)./sum(W,1);
end

%{
estimation of standard deviation with weights
using Rousseeuw and Croux method
for each row, calculate std of columns
%}
function s = std_robust_w(X,W)
X = double(X);
W = double(W);
triul = @(t)reshape(t(~diag(ones(1,size(t, 1)))), size(t)-[1 0]);
s = zeros(1,size(X,2));
for j = 1:size(X,2),
%     med = median_w(X(:,j),W(:,j));
%     v = sqrt(var(dat(:,j)));
%     v = 1.4826*median(pdist(dat(:,j),'cityblock'));
    v = squareform( pdist(X(:,j),'cityblock') );
    w = W(:,j)*W(:,j)';
    v = triul(v);
    w = triul(w);
%     v_col_med = median_w(v,w);
    v_col_med = median(v);
    col_w = zeros(1,size(w,2));
    for k = 1:size(w,2),
        col_w(k) = mean(w(v(:,k) == v_col_med(k),k));
    end
%     s(j) = 1.1926*median_w(v_col_med',col_w');
    s(j) = 1.1926*median(v_col_med');
end
end

%{
estimation of standard deviation with weights
for each row, calculate std of columns

if y = X*a', then calculate std of rows of y
%}
function s = std_w(X,W,a)
X = double(X);
W = double(W);

is_a = 1;
if ~exist('a', 'var')
    is_a = 0;
end

if (is_a),
    a = double(a);
	%example, s = std_w(X,ones(size(X)),pc_w(:,i_pc)')
    if size(X,1) == size(a, 2),
        X = X';
        W = W';
        a = a';        
%         mode = 4;
%         rho = pdet_util_c(X,W,mode,a);

        %exact center       
        [n,p] = size(X);
        [~,m] = size(a);
        
        tic
        s = zeros(1,n);
        k = 1;
        for i = 1:n,
            xw = zeros(1,m);
%             yw = zeros(n-i,m);
%             xy = zeros(n-i,1);
            xx = zeros(1,1);
%             yy = zeros(n-i,1);

            u_x = zeros(1,1);
%             u_y = zeros(n-i,1);

            w_sum = zeros(1,1);

            %speed up indexing
            x_i = X(i,:);
%             y_i = X((i+1):n,:);
            w_x_i = W(i,:);
%             w_y_i = W((i+1):n,:);

            for q = 1:p,    %matlab jit
                w_j = w_x_i(q);
%                 w_j = w_x_i(q)*w_y_i(:,q);
                w_sum = w_sum + w_j;
                x = x_i(q).*w_j;
%                 y = y_i(:,q).*w_j;        
                for kk = 1:m,
                    xw(:,kk) = xw(:,kk) + x.*a(q,kk);
%                     yw(:,kk) = yw(:,kk) + y.*a(q,kk);
                end
            end

            for kk = 1:m,
                u_x = u_x + xw(:,kk);
%                 u_y = u_y + yw(:,kk);
            end
            w_sum = w_sum/p;
            u_x = u_x./w_sum/m;
%             u_y = u_y./w_sum/m;
%             u_x = u_x/m;
%             u_y = u_y/m;
            for kk = 1:m,
                x = xw(:,kk)./w_sum - u_x;
%                 y = yw(:,kk)./w_sum - u_y;
%                 xy = xy + x.*y;
                xx = xx + x.*x;
%                 yy = yy + y.*y;
%                 xy = xy + xw(:,kk).*yw(:,kk);
%                 xx = xx + xw(:,kk).*xw(:,kk);
%                 yy = yy + yw(:,kk).*yw(:,kk);                
            end
%             xy = xy - u_x.*u_y*m;
%             xx = xx - u_x.*u_x*m;
%             yy = yy - u_y.*u_y*m;
            
            s(k) = sqrt(xx/(m-1));
            k = k + 1;
            if rem(i,100) == 0,
                disp(sprintf('iteration %d',i));
                toc
            end
        end
    else
        %approximate center
        %exact center would use pairwise w_j
        %example, s = std_w(X',ones(size(X')),pc_w(:,i_pc)');
        c = cov_w(X, W, a);
        s = sqrt(diag(c));
    end
else
    X = zeromean_w(X,W);
    s = sqrt(sum(X.^2.*W,1)./sum(W,1)*size(X,1)/(size(X,1)-1));
end

end

%{
pairwise mean of elements of a vector
%mode: {'geomean1','geomean0'}
%}
function s = pmean_w(X,W,mode)
X = double(X);
W = double(W);
n = length(X);
if strcmp(mode,'geomean0'),
    X = log(X+realmin);
elseif strcmp(mode,'geomean1'),
    X = log(X+1);
end

s = zeros(1,n*(n-1)./2);
k = 1;
for i = 1:n,
    x_i = X(i);
    y_i = X((i+1):n);
    w_x_i = W(i);
    w_y_i = W((i+1):n);
    
    s(k:(k+n-i-1)) = (x_i.*w_x_i + y_i.*w_y_i)./(w_x_i + w_y_i);
    k = k + (n-i);
end

if strcmp(mode,'geomean0'),
    s = exp(s);
elseif strcmp(mode,'geomean1'),
    s = exp(s)-1;
end

end

%{
prepare Ec and Er for isa2
X: row: observations, column: variables
%}
function [Ec,Er] = prepEcEr_isa2_w(X,W)
X = double(X);
W = double(W);

var_adj_Ec = 0.7;
var_adj_Er = 0.7;

[n,p] = size(X);

Ec = double(X);
Er = double(X');

i_z = std(Ec,[],1) == 0;
Ec(:,i_z) = [];
Er(i_z,:) = [];
W(:,i_z) = [];

%adjust variance
Ec = bsxfun(@minus, Ec, mean_w(Ec,W)); Ec = bsxfun(@rdivide, Ec, std_w(Ec,W).^(var_adj_Ec)); %genes are coordinates, cells are observations
if false,
    Er = Ec';
end
Er = bsxfun(@minus, Er, mean_w(Er,W')); Er = bsxfun(@rdivide, Er, std_w(Er,W').^(var_adj_Er)); %cells are coordinates, genes are observations    

Ec = Ec.*W;
Er = Er.*W';
end

%{
%calculate spatial cross correlation with w
%between genes in g1 and genes in g2
%ydata comes from tsne, should have local dense regions
%option = {'moranI',@tpdf,1,g1,g2}
function p = spatialcrosscorr_w(X,W,ydata,option)
X = double(X);
W = double(W);

mode = option{1};

if strcmp(mode,'moranI'),
    qpdf = option{2};
    kk = option{3};  %select for gaussian or randomization
    g1 = option{4};
    g2 = option{5};
    
    lg1 = length(g1);
    lg2 = length(g2);
    g12 = intersect(g1,g2);
     
    %use Moran's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0; %how much corr versus spatiall corr
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    EI = -1/(n-1);
%     S0 = 2 * sum(Q(:));
%     S1 = 2 * sum(Q(:).^2);
%     S2 = 4 * sum(sum(Q,1).^2);
%     A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
%     B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
%     B = S1 - 2*n*S2 + 6*S0^2;
%     C = (n-1)*(n-2)*(n-3)*S0^2;
    
    p = zeros(lg1,lg2) + EI;
%     V = ones(lg1,lg2);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    tic
    for ij = 1:lg1,
        for ijj = 1:lg2,
            j = g1(ij);
            jj = g2(ijj);

            if j == jj,
                continue
            end
            if ismember(j,g12) && ismember(jj,g12),
                if p(ij,ijj) ~= EI,
                    p(ijj,ij) = p(ij,ijj);
%                     V(ijj,ij) = V(ij,ijj);
                    continue
                elseif p(ijj,ij) ~= EI,
                    p(ij,ijj) = p(ijj,ij);
%                     V(ij,ijj) = V(ijj,ij);
                    continue
                end
            end

            x = double(X(:,j)');
            wx = double(W(:,j)');
            x = x - mean_w(x',wx');

            y = double(X(:,jj)');
            wy = double(W(:,jj)');
            y = y - mean_w(y',wy');

            if ~sum(x) || ~sum(y),
                continue
            end

%             w = wx.*wy;
%             D = ( sum((x2.*y2) .* w) / sum(w) )/...
%                 (( sum(x2.*w)/sum(w) ) * ( sum(y2.*w)/sum(w) ));

%             EI2 = (A - D*B)/C;
%             V(ij,ijj) = EI2 - EI^2;

            %Q off diagonal elements
            xy_sum = 0;
            xx_sum = 0;
            yy_sum = 0;
            s0 = 0;
            k = 1;
            for i = 1:n,
                wq = (wx(i).*wy(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((x(i)*y(i+1:n)) .* wq);
                xx_sum = xx_sum + sum((x(i)*x(i+1:n)) .* wq);
                yy_sum = yy_sum + sum((y(i)*y(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end
            k = 1;
            for i = 1:n,
                wq = (wy(i).*wx(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((y(i)*x(i+1:n)) .* wq);
                xx_sum = xx_sum + sum((x(i)*x(i+1:n)) .* wq);
                yy_sum = yy_sum + sum((y(i)*y(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end

            p(ij,ijj) = (xy_sum) / ...
                sqrt( abs(xx_sum * yy_sum) );
        end
        if ~rem(ij,10),
            fprintf('ij = %d, ',ij);
            toc
        end
    end
    
    %randomization variance sometimes gives erroneous result
%     Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
%     if kk,
%         V = Vg;
%     else
%         fprintf('V<0: %d, V: %d',sum(V<0),length(V))
%         V(V<Vg) = Vg;
%     end
%     pz = -abs((p - EI)./sqrt(V));
%     pz = 2*normcdf(pz);
%     p = pz;
else,
    p = 1;
end
end
%}

%{
calculate pairwise spatial cross correlation with weights between genes in list1 and in list2
between genes in g1 and genes in g2

ydata: tsne 2d embedding
option.mode = 'moranI','moranI_th','gearyC'
option.pdf = @tpdf, @normpdf
option.g1 = list1 of gene index
option.g2 = list2 of gene index
%}
function p = spatialcorr_w(X,W,ydata,option)
X = double(X);
W = double(W);

mode = option.mode;

if strcmp(mode,'moranI'),
    qpdf = option.pdf;
%     kk = option{3};  %select for gaussian or randomization
    g1 = option.g1;
    g2 = option.g2;
    
    lg1 = length(g1);
    lg2 = length(g2);
    g12 = intersect(g1,g2);
     
    %use Moran's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) );

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0; %how much corr versus spatiall corr
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    EI = -1/(n-1);
%     S0 = 2 * sum(Q(:));
%     S1 = 2 * sum(Q(:).^2);
%     S2 = 4 * sum(sum(Q,1).^2);
%     A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
%     B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
%     B = S1 - 2*n*S2 + 6*S0^2;
%     C = (n-1)*(n-2)*(n-3)*S0^2;
    
    p = zeros(lg1,lg2) + EI;
%     V = ones(lg1,lg2);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    tic
    for ij = 1:lg1,
        for ijj = 1:lg2,
            j = g1(ij);
            jj = g2(ijj);

            if j == jj,
                continue
            end
            if ismember(j,g12) && ismember(jj,g12),
                if p(ij,ijj) ~= EI,
                    p((g1 == jj), (g2 == j)) = p(ij,ijj);
%                     V(ijj,ij) = V(ij,ijj);
                    continue
                elseif p((g1 == jj), (g2 == j)) ~= EI,
                    p(ij,ijj) = p((g1 == jj), (g2 == j));
%                     V(ij,ijj) = V(ijj,ij);
                    continue
                end
            end

            x = double(X(:,j)');
            wx = double(W(:,j)');
            x = x - mean_w(x',wx');
            x2 = x.^2;

            y = double(X(:,jj)');
            wy = double(W(:,jj)');
            y = y - mean_w(y',wy');
            y2 = y.^2;

            if sum(x(x>0)) < 1e-10 || sum(y(y>0)) < 1e-10,
                continue
            end

            w = wx.*wy;
%             D = ( sum((x2.*y2) .* w) / sum(w) )/...
%                 (( sum(x2.*w)/sum(w) ) * ( sum(y2.*w)/sum(w) ));

%             EI2 = (A - D*B)/C;
%             V(ij,ijj) = EI2 - EI^2;

            %Q off diagonal elements
            xy_sum = 0;
            s0 = 0;
            k = 1;
            for i = 1:n,
                wq = (wx(i).*wy(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((x(i)*y(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end
            k = 1;
            for i = 1:n,
                wq = (wy(i).*wx(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((y(i)*x(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end

            p(ij,ijj) = (xy_sum/s0) / ...
                sqrt( sum(x2.*w)/sum(w) * sum(y2.*w)/sum(w) );
        end
        if ~rem(ij,10),
            fprintf('ij = %d, ',ij);
            toc
        end
    end
    
    %transform moran I to be symmetric and within [-1,1]
    %Maruyama 2015
%     if length(option)>= 5,
%         if option{5},
            %construct H, Helmert matrix
            H = zeros(n,n-1);
            for i = 1:n-1,
                H(1:i+1,i) = [ones(1,i),-i]'/sqrt(i*(i+1));
            end
            Q = squareform(Q);
            W_tilde = H'*Q*H/(sum(Q(:))/n);
            [~,D] = eig(W_tilde);
            lambda_1 = min(diag(D));
            lambda_n = max(diag(D));
            disp(sprintf('min lambda : %.2f',lambda_1));
            disp(sprintf('max lambda : %.2f',lambda_n));
            pp = (n-1)*p+1;
            pp(pp<0) = pp(pp<0)/abs( (n-1)*lambda_1 + 1 );
            pp(pp>=0) = pp(pp>=0)/abs( (n-1)*lambda_n + 1 );
            p = {p,pp};
%         end
%     end
            
        
    %randomization variance sometimes gives erroneous result
%     Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
%     if kk,
%         V = Vg;
%     else
%         fprintf('V<0: %d, V: %d',sum(V<0),length(V))
%         V(V<Vg) = Vg;
%     end
%     pz = -abs((p - EI)./sqrt(V));
%     pz = 2*normcdf(pz);
%     p = pz;

elseif strcmp(mode,'moranI_th'),
    qpdf = option.pdf;
%     kk = option{3};  %select for gaussian or randomization
    g1 = option.g1;
    g2 = option.g2;
    
    lg1 = length(g1);
    lg2 = length(g2);
    g12 = intersect(g1,g2);
     
    %use Moran's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) );

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0; %how much corr versus spatiall corr
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    EI = -1/(n-1);
%     S0 = 2 * sum(Q(:));
%     S1 = 2 * sum(Q(:).^2);
%     S2 = 4 * sum(sum(Q,1).^2);
%     A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
%     B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
%     B = S1 - 2*n*S2 + 6*S0^2;
%     C = (n-1)*(n-2)*(n-3)*S0^2;
    
    p = zeros(lg1,lg2) + EI;
%     V = ones(lg1,lg2);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    X = double(X>1.1);
    
    tic
    for ij = 1:lg1,
        for ijj = 1:lg2,
            j = g1(ij);
            jj = g2(ijj);

            if j == jj,
                continue
            end
            if ismember(j,g12) && ismember(jj,g12),
                if p(ij,ijj) ~= EI,
                    p((g1 == jj), (g2 == j)) = p(ij,ijj);
%                     V(ijj,ij) = V(ij,ijj);
                    continue
                elseif p((g1 == jj), (g2 == j)) ~= EI,
                    p(ij,ijj) = p((g1 == jj), (g2 == j));
%                     V(ij,ijj) = V(ijj,ij);
                    continue
                end
            end

            x = double(X(:,j)');
            wx = double(W(:,j)');
            x = x - mean_w(x',wx');
            x2 = x.^2;

            y = double(X(:,jj)');
            wy = double(W(:,jj)');
            y = y - mean_w(y',wy');
            y2 = y.^2;

            if sum(x(x>0)) < 1e-10 || sum(y(y>0)) < 1e-10,
                continue
            end

            w = wx.*wy;
%             D = ( sum((x2.*y2) .* w) / sum(w) )/...
%                 (( sum(x2.*w)/sum(w) ) * ( sum(y2.*w)/sum(w) ));

%             EI2 = (A - D*B)/C;
%             V(ij,ijj) = EI2 - EI^2;

            %Q off diagonal elements
            xy_sum = 0;
            s0 = 0;
            k = 1;
            for i = 1:n,
                wq = (wx(i).*wy(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((x(i)*y(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end
            k = 1;
            for i = 1:n,
                wq = (wy(i).*wx(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((y(i)*x(i+1:n)) .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end

%             p(ij,ijj) = (xy_sum/s0) / ...
%                 sqrt( sum(x2.*w)/sum(w) * sum(y2.*w)/sum(w) );
            p(ij,ijj) = (xy_sum/s0);            
        end
        if ~rem(ij,10),
            fprintf('ij = %d, ',ij);
            toc
        end
    end
    
    %transform moran I to be symmetric and within [-1,1]
    %Maruyama 2015
%     if length(option)>= 5,
%         if option{5},
            %construct H, Helmert matrix
            H = zeros(n,n-1);
            for i = 1:n-1,
                H(1:i+1,i) = [ones(1,i),-i]'/sqrt(i*(i+1));
            end
            Q = squareform(Q);
            W_tilde = H'*Q*H/(sum(Q(:))/n);
            [~,D] = eig(W_tilde);
            lambda_1 = min(diag(D));
            lambda_n = max(diag(D));
            disp(sprintf('min lambda : %.2f',lambda_1));
            disp(sprintf('max lambda : %.2f',lambda_n));
            pp = (n-1)*p+1;
            pp(pp<0) = pp(pp<0)/abs( (n-1)*lambda_1 + 1 );
            pp(pp>=0) = pp(pp>=0)/abs( (n-1)*lambda_n + 1 );
            p = {p,pp};
%         end
%     end
            
        
    %randomization variance sometimes gives erroneous result
%     Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
%     if kk,
%         V = Vg;
%     else
%         fprintf('V<0: %d, V: %d',sum(V<0),length(V))
%         V(V<Vg) = Vg;
%     end
%     pz = -abs((p - EI)./sqrt(V));
%     pz = 2*normcdf(pz);
%     p = pz;

elseif strcmp(mode, 'gearyC')   %do not use
    qpdf = option.pdf;
%     kk = option{3};  %select for gaussian or randomization
    g1 = option.g1;
    g2 = option.g2;
    
    lg1 = length(g1);
    lg2 = length(g2);
    g12 = intersect(g1,g2);
     
    %use Geary's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0;   %not counting self
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);
    
    EC = 1;
%     S0 = 2 * sum(Q(:));
%     S1 = 2 * sum(Q(:).^2);
%     S2 = 4 * sum(sum(Q,1).^2);
%     C = n*(n-2)*(n-3)*S0^2;

    p = zeros(lg1,lg2) + EC;
%     V = ones(lg1,lg2);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    tic
    for ij = 1:lg1,
        for ijj = 1:lg2,
            j = g1(ij);
            jj = g2(ijj);

            if j == jj,
                continue
            end
            if ismember(j,g12) && ismember(jj,g12),
                if p(ij,ijj) ~= EC,
                    p(ijj,ij) = p(ij,ijj);
                    continue
                elseif p(ijj,ij) ~= EC,
                    p(ij,ijj) = p(ijj,ij);
                    continue
                end
            end

            x = double(X(:,j)');
            wx = double(W(:,j)');
            zx2 = (x - mean_w(x',wx')).^2;

            y = double(X(:,jj)');
            wy = double(W(:,jj)');
            zy2 = (y - mean_w(y',wy')).^2;

            if ~sum(x) || ~sum(y),
                continue
            end
            
            w = wx.*wy;
%             D = ( sum((zx2.*zy2) .* w) / sum(w) )/...
%                 (( sum(zx2.*w)/sum(w) ) * ( sum(zy2.*w)/sum(w) ));
            
%             A = (n-1)*S1*(n^2-3*n+3-(n-1)*D) +...
%                 S0^2*(n^-3-(n-1)^2*D);
%             B = 1/4*(n-1)*S2*(n^2+3*n-6-(n^2-n+2)*D);
%             V(ij,ijj) = (A-B)/C;

            %Q off diagonal elements
            xy_sum = 0;
            s0 = 0;
            k = 1;
            for i = 1:n,
                wq = (wx(i).*wy(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((x(i)-y(i+1:n)).^2 .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end
            k = 1;
            for i = 1:n,
                wq = (wy(i).*wx(i+1:n)) .* Q(k:(k+n-i-1));
                xy_sum = xy_sum + sum((y(i)-x(i+1:n)).^2 .* wq);
                s0 = s0 + sum(wq);
                k = k + (n-i);
            end

            %Q diagonal elements should be zeros
            p(ij,ijj) = 1/2 * (xy_sum/s0) /...
                sqrt( sum( zx2.*w )/(sum(w)-1) *...
                sum( zy2.*w )/(sum(w)-1) );
        end
        if ~rem(ij,10),
            fprintf('ij = %d, ',ij);
            toc
        end
    end
    
%     Vg = 1/( 2*(n+1)*S0.^2 ) * ((2*S1 + S2)*(n-1) - 4*S0.^2);    %gaussian    
%     if kk,
%         V = Vg;
%     else
%         fprintf('V<0: %d, V: %d\n',sum(V<0),length(V))
%         V(V<Vg) = Vg;
%     end

%     pz = -abs((p - EC)./sqrt(V));
%     pz = 2*normcdf(pz);
%     p = pz;
else
    p = 1;
end

end

%{
two-sample likelihood ratio test
using combined binomial and normal model

for each row, test across columns
rows of x corresponds to observations
%}
function [H,p] = lratiotest2_binorm_w(X,Y,W_X,W_Y)
X = double(X); Y = double(Y);
W_X = double(W_X); W_Y = double(W_Y);
samples = {X,Y,[X;Y]};
samples_w = {W_X,W_Y,[W_X;W_Y]};

logl = cell(1,length(samples));
p_bi = cell(1,2);
p_norm = cell(1,2);
for i = 1:length(samples),
    %fit binomial distribution
    X = samples{i};
    W_X = samples_w{i};
    i_nz = X > 0;
    i_nz_low = X>0 & X<2.1;
    i_nz_high = X>2.1;
    n = size(X,1);
    k = (sum(i_nz,1) + sum((1-W_X).*(~i_nz),1));
    k = min(max(k,1),n);
    if i < 3,
        p = k/n;
        l_bi = binopdf(ceil(k),n,p);
        p_bi{i} = p;
    else
        p = p_bi{1};
        l_bi_1 = binopdf(ceil(k),n,p);
        p = p_bi{2};
        l_bi_2 = binopdf(ceil(k),n,p);
        l_bi = max([l_bi_1;l_bi_2]);
    end
    
    %fit normal distribution
    l_norm = zeros(1,size(X,2));
    ux = zeros(1,size(X,2));
    sx = zeros(1,size(X,2));
    if i < 3,
        for j = 1:size(X,2),
            if sum(X(i_nz_high(:,j),j)) == 0 || length(X(i_nz_high(:,j),j)) < 3,
                ux_j = 0;
                sx_j = 0;
                l_norm(j) = 1;
            else
                ux_j = mean(X(i_nz_high(:,j),j));
                sx_j = std(X(i_nz_high(:,j),j));
                l_norm(j) = prod(normpdf(X(i_nz_high(:,j),j),ux_j,sx_j));
            end
            ux(j) = ux_j;
            sx(j) = sx_j;
        end
        p_norm{i} = {ux,sx};
    else
        for j = 1:size(X,2),
            if sum(X(i_nz_high(:,j),j)) == 0 || length(X(i_nz_high(:,j),j)) < 3,
                l_norm(j) = 1;
            else
                ux_j = p_norm{1}{1}(j);
                sx_j = p_norm{1}{2}(j);
                l_norm_1 = prod(normpdf(X(i_nz_high(:,j),j),ux_j,sx_j));
                ux_j = p_norm{2}{1}(j);
                sx_j = p_norm{2}{2}(j);
                l_norm_2 = prod(normpdf(X(i_nz_high(:,j),j),ux_j,sx_j));
                l_norm(j) = max(l_norm_1, l_norm_2);
            end
        end
    end    
    logl{i} = log(l_bi.*l_norm);
end

dof = 3;
alpha = 0.05;
uLL = logl{1} + logl{2};
rLL = logl{3};
stat = 2*(uLL-rLL);
p = 1-chi2cdf(stat,dof);
H = (p <= alpha);

end

%{
two-sample binomial and normal test

for each row, test across columns
rows of x corresponds to observations
%}
function [H,p] = zttest2_w(X,Y,W_X,W_Y)
X = double(X); Y = double(Y);
W_X = double(W_X); W_Y = double(W_Y);
samples = {X,Y,[X;Y]};
samples_w = {W_X,W_Y,[W_X;W_Y]};

p_bi = cell(1,2);
for i = 1:length(samples),
    %fit binomial distribution
    X = samples{i};
    W_X = samples_w{i};
    i_nz = X > 0;
%     i_nz_low = X>0 & X<2.1;
%     i_nz_high = X>2.1;
    n = size(X,1);
    k = (sum(i_nz,1) + sum((1-W_X).*(~i_nz),1));
    p_bi{i} = {min(max(k/n,1/n),(n-1)/n),n};
    samples_w{i} = double(i_nz);
end

% p1 = min(max(p_bi{1}/n,1/n),(n-1)/n);
% p2 = min(max(p_bi{2}/n,1/n),(n-1)/n);
% k1 = min(ceil(p_bi{1}),n);
% k2 = min(ceil(p_bi{2}),n);
% % l_bi = min([binopdf(k1,n,p2); binopdf(k2,n,p1)]);
% l_bi = min([binopdf(k1,n,p2); binopdf(k2,n,p1)]);
p1 = p_bi{1}{1};
p2 = p_bi{2}{1};
p = p_bi{3}{1};
n1 = p_bi{1}{2};
n2 = p_bi{2}{2};
z = (p1 - p2) ./ sqrt( p.*(1-p) * (1/n1 + 1/n2));
l_bi = 2*normcdf( -abs(z) );

[~,l_norm] = ttest2_w(samples{1},samples{2},samples_w{1},samples_w{2});
l_norm_z = sum(samples_w{1},1) < 3 | sum(samples_w{2},1) < 3;
l_norm(l_norm_z) = 1;
p = min([l_bi;l_norm]);
H = p<0.05;

end
