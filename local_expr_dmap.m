%{
estimate how rare it is to find a gene highly expressed in cells 
positioned in a neighborhood in tSNE embedding

Yinqing Li
yinqing@csail.mit.edu
2015

inputs:
X: rows are observations (cells), columns are variables (genes)
W: weights
ydata: tsne 2d embedding

outputs:
p: p value for variables

option:

option.mode = 'hg','rand','ostat','moranI','moranI_th','gearyC','aple'

options for different statistical methods

hypergeometric test
assume independence between pairs of cells
option.pdf: probability density function for transforming pairwise distances to similarities
option.k: expected number of close neighbors
local_expr_dmap(X',W',ydata,struct('mode','hg','pdf',@tpdf,'k',6));

randomization test
use randomization to find the background distribution
option.pdf: probability density function for transforming pairwise distances to similarities
option.k: expected number of minimal close neighbors
option.n_rand: number of random sampling
local_expr_dmap(X',W',ydata,struct('mode','rand','pdf',@tpdf,'k',6,'n_rand',1e5));

order statistics
use manhattan distance
option.gp_th: significance level for removing outliers tested using normal distribution
option.th_qtle: whether to use 0.1 quantile as the threshold to binarize expression level
local_expr_dmap(X',W',ydata,struct('mode','ostat','gp_th',0.02,'th_qtle',1.1));

Moran's I
option.pdf: probability density function for transforming pairwise distances to similarities
option.prox: whether to use approximated background distribution in calculation of p value
local_expr_dmap(X',W',ydata,struct('mode','moranI','pdf',@(x) tpdf(x,1),'prox',1));

Moran's I with a threshold
option.pdf: probability density function for transforming pairwise distances to similarities
option.prox: whether to use approximated background distribution in calculation of p value
local_expr_dmap(X',W',ydata,struct('mode','moranI_th','pdf',@(x) normpdf(x,0,2),'prox',1));

Geary's C
option.pdf: probability density function for transforming pairwise distances to similarities
option.prox: whether to use approximated background distribution in calculation of p value
local_expr_dmap(X',W',ydata,struct('mode','gearyC','pdf',@(x) tpdf(x,1),'prox',1));

APLE
option.pdf: probability density function for transforming pairwise distances to similarities
option.prox: whether to use approximated background distribution in calculation of p value
local_expr_dmap(X',W',ydata,struct('mode','aple','pdf',@(x) tpdf(x,1),'prox',1));

Copyright (c) 2015, Yinqing Li
%}

function p = local_expr_dmap(X,W,ydata,option)

fh = func_util_bisne;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

fh = func_util_pdet;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

mode = option.mode;

X = double(X);
W = double(W);

if strcmp(mode, 'hg'),  %hypergeometric, works for small n
    
    qpdf = option.pdf;
    k = option.k;
    
    %estimate t distribution and probability for each sample
    %t distribution with 1 degree of freedom
    Q = qpdf( (pdist(ydata,'euclidean')) , 1);

    %normalize Q
    Q = max(Q / sum(Q(:)), realmin);
    [n,m] = size(X);

    if isfield(option,'n_sc'),
        n_sc = option.n_sc;
        if ~n_sc || n_sc>n,
            n_sc = n;
        end
    else
        n_sc = n;
    end
    a = (n_sc/n)^2;
    
    %use hypergeometric test for pairwise distance and expression
    %note that pairwise metrics are not completely independent
    %thresholding expression and distance
    th_Q = quantile(Q(:),1-k/(n-1));
    idx_Q = Q>th_Q;
    K = sum(idx_Q);
    M = length(Q);
    
    K = floor(K*a);
    M = floor(M*a);
            
    p = ones(1,m);
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        if sum(x(x>0)) < 1e-10,
            continue
        end

    %     %thresholding method 1
    %     x = double(x>1);
    %     P = pmean_w(x,w,'geomean0');
    %     th_P = 0.5;

        %thresholding method 2
        P = pmean_w(x,w,'geomean0');
        th_P = quantile(x(x>0),0.1);  

        idx_P = P>th_P;
        x = sum(idx_P & idx_Q);
        N = sum(idx_P);
        
        x = floor(x*a);
        N = floor(N*a);
        if ~x || ~N,
            continue
        end
        
        p(j) = hygecdf(x,M,K,N,'upper');

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    
elseif strcmp(mode, 'rand'),    %randomization test

    qpdf = option.pdf;
    k = option.k;
    if isfield(option,'n_rand'),
        n_rand = option.n_rand;
    else,
        n_rand = 1e6;
    end
    
    %estimate t distribution and probability for each sample
    %t distribution with 1 degree of freedom
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);
    [n,m] = size(X);
    
    %normalize Q
    Q(1:n+1:end) = 0;
    Q = max(Q / sum(Q(:)), realmin);
    [n,m] = size(X);

    %use hypergeometric test for pairwise distance and expression
    %note that pairwise metrics are not completely independent
    %thresholding expression and distance
    th_Q = quantile(Q(:),1-k/(n-1));
    idx_Q = Q>th_Q;

    K = sum(idx_Q);
    M = length(Q);
    
    idx_Q = squareform(idx_Q);
    
    sum_k = zeros(n,n_rand);
    sprintf('estimated time is %.02fs',5e-7*n*n_rand)
    tic
    for jj = 1:n_rand,
        ii = randperm(n);
        P = tril(idx_Q(ii,ii));
        sum_k(:,jj) = sum(P,2);
    end
    sum_k = cumsum(sum_k,1);
    %keep on the top 1000 randomization
    n_rand_s = 1e3;
    sum_k_s = zeros(n,n_rand_s);
    for jj = 1:n,
        Y = sort(sum_k(jj,:),'descend');
        sum_k_s(jj,:) =Y(1:n_rand_s);
    end
    sum_k = sum_k_s;
    clear sum_k_s n_rand_s
	toc    
    
    idx_Q = squareform(idx_Q);
    
    p = ones(1,m);
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        if sum(x(x>0)) < 1e-10,
            continue
        end
        
        %thresholding method 2
        P = pmean_w(x,w,'geomean0');
%         P = squareform(pmean_w(x,w,'geomean0'));
%         P(1:n+1:end) = 0;
        th_P = quantile(x(x>0),0.1);  

        idx_P = P>th_P;
        k = sum(sum(idx_P & idx_Q));
        N = sum(x>th_P);
        
        if ~N || ~k,
            continue
        end
        
        p(j) = sum(sum_k(N, sum_k(N,:)>=k))/n_rand;
        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end        
    end
    
    %{
    n_rand = 1e4;

    p = ones(1,m);
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        if ~sum(x>0)
            continue
        end

    %     %thresholding method 1
    %     x = double(x>1);
    %     P = pmean_w(x,w,'geomean0');
    %     th_P = 0.5;

        %thresholding method 2
        P = squareform(pmean_w(x,w,'geomean0'));
        P(1:n+1:end) = 0;
        th_P = quantile(x(x>0),0.1);  

        idx_P = P>th_P;
        x = sum(sum(idx_P & idx_Q));
        N = sum(idx_P(:));

        if ~x || ~N,
            continue
        end

        %randomization test
        c_x_perm = 0;
        x_perm_sum = 0;
        x_perm2_sum = 0;
        for jj = 1:n_rand,
            ii = randperm(n);
            x_perm = sum(sum(idx_P(ii,ii) & idx_Q));
            x_perm_sum = x_perm_sum + x_perm;
            x_perm2_sum = x_perm2_sum + x_perm^2;
            c_x_perm = c_x_perm + x_perm > x;
        end
        
        %approximation
        if if_gaussian
            mean_x_perm = x_perm_sum/n_rand;
            std_x_perm = sqrt((x_perm2_sum - mean_x_perm^2*n_rand)/(n_rand-1));
            z = (x - mean_x_perm)/std_x_perm;
            p(j) = 2*normcdf(-abs(z),0,1);
        else
            p(j) = c_x_perm/n_rand;
        end
        
        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    %}

elseif strcmp(mode, 'ostat'),  %order statistics
    
%     fcom = option.fcom;   %method for combining p values
    gp_th = option.gp_th;  %threshold for outlier
    if_quantile = option.th_qtle;   %use 0.1 quantile as threshold
    
    ydata = bsxfun(@minus,ydata,min(ydata,[],1) );
    ydata = bsxfun(@rdivide,ydata,max(ydata,[],1));
    
    [~,m] = size(X);
    
    pp = ones(m,size(ydata,2));
    tic
    for j = 1:m,
        x = double(X(:,j)');
%         w = double(W(:,j)');
        if sum(x(x>0)) < 1e-10,
            continue
        end

    %     %thresholding method 1
    %     x = double(x>1);
    %     P = pmean_w(x,w,'geomean0');
    %     th_P = 0.5;

        %thresholding method 2
        th_x = quantile(x(x>0),0.1);
        i = x>th_x;
        if ~if_quantile,
            i = x>0;
        end
        n = sum(i);
        
        if n<2,
            continue
        end
        
        pp(j,:) = ostat_w(ydata(i,:), gp_th);

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    p = pp;
    
elseif strcmp(mode, 'moranI'),
    
    qpdf = option.pdf;
    kk = option.prox;  %select for gaussian or randomization
     
    %use Moran's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) );

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0;   %not counting self
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    EI = -1/(n-1);
    S0 = 2 * sum(Q(:));
    S1 = 2 * sum(Q(:).^2);
    S2 = 4 * sum(sum(Q,1).^2);
    A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
    B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
    C = (n-1)*(n-2)*(n-3)*S0^2;
    
    pp = zeros(1,m) + EI;
    V = ones(1,m);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        x = x - mean_w(x',w');
        
        if sum(x(x>0)) < 1e-10,
            continue
        end
        
        D = ( sum((x.^4).*w) / sum(w) )/...
            ( sum((x.^2).*w) / sum(w) ).^2;
        EI2 = (A - D*B)/C;
        V(j) = EI2 - EI^2;
        
        %Q off diagonal elements
        x_sum = 0;
        s0 = 0;
        k = 1;
        for i = 1:n,
            wq = (w(i).*w(i+1:n)) .* Q(k:(k+n-i-1));
            x_sum = x_sum + sum((x(i)*x(i+1:n)) .* wq);
            s0 = s0 + sum(wq);
            k = k + (n-i);
        end
        x_sum = x_sum*2;
        s0 = s0*2;
        
        %Q diagonal elements
        pp(j) = (x_sum/s0) / (sum(x.^2.*w)/sum(w));

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    
    %randomization variance sometimes gives erroneous result
    Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
    if kk,
        V = Vg;
    else
        fprintf('V<0: %d, V<Vg: %d, V: %d\n',sum(V<0),sum(V<Vg),length(V))
        V(V<Vg) = Vg;
    end
    pz = - ((pp - EI)./sqrt(V));    %large pp gets small p value
    pz = 2*normcdf(pz);
    p = pz;

elseif strcmp(mode, 'moranI_th'),
    
    qpdf = option.pdf;
    kk = option.prox;  %select for gaussian or randomization
     
    %use Moran's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) );

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0;   %not counting self
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    EI = -1/(n-1);
    S0 = 2 * sum(Q(:));
    S1 = 2 * sum(Q(:).^2);
    S2 = 4 * sum(sum(Q,1).^2);
    A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
    B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
    C = (n-1)*(n-2)*(n-3)*S0^2;
    
    pp = zeros(1,m) + EI;
    V = ones(1,m);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    X = double(X>1.1);
    
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        x = x - mean_w(x',w');
        
        if sum(x(x>0)) < 1e-10,
            continue
        end
        
        D = ( sum((x.^4).*w) / sum(w) )/...
            ( sum((x.^2).*w) / sum(w) ).^2;
        EI2 = (A - D*B)/C;
        V(j) = EI2 - EI^2;
        
        %Q off diagonal elements
        x_sum = 0;
        s0 = 0;
        k = 1;
        for i = 1:n,
            wq = (w(i).*w(i+1:n)) .* Q(k:(k+n-i-1));
            x_sum = x_sum + sum((x(i)*x(i+1:n)) .* wq);
            s0 = s0 + sum(wq);
            k = k + (n-i);
        end
        x_sum = x_sum*2;
        s0 = s0*2;
        
        %Q diagonal elements
        pp(j) = (x_sum/s0);
%         pp(j) = (x_sum/s0) / (sum(x.^2.*w)/sum(w));

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    
    %randomization variance sometimes gives erroneous result
    Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
    if kk,
        V = Vg;
    else
        fprintf('V<0: %d, V<Vg: %d, V: %d\n',sum(V<0),sum(V<Vg),length(V))
        V(V<Vg) = Vg;
    end
    pz = - ((pp - EI)./sqrt(V));    %large pp gets small p value
    pz = 2*normcdf(pz);
    p = pz;
    
elseif strcmp(mode, 'gearyC'),
    
    qpdf = option.pdf;
    kk = option.prox;
     
    %use Geary's I
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0;   %not counting self
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);
    
    EC = 1;
    S0 = 2 * sum(Q(:));
    S1 = 2 * sum(Q(:).^2);
    S2 = 4 * sum(sum(Q,1).^2);
    C = n*(n-2)*(n-3)*S0^2;
    
    pp = zeros(1,m) + EC;
    V = ones(1,m);
    
    Q(1:n+1:end) = 0;
    Q = squareform(Q);
    
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        ux = mean_w(x',w');
        
        if sum(x>10) < 1e-10,
            continue
        end
        
        D = ( sum((x-ux).^4.*w) / sum(w) )/...
            ( sum((x-ux).^2.*w) / sum(w) ).^2;
        A = (n-1)*S1*(n^2-3*n+3-(n-1)*D) +...
            S0^2*(n^-3-(n-1)^2*D);
        B = 1/4*(n-1)*S2*(n^2+3*n-6-(n^2-n+2)*D);
        V(j) = (A-B)/C;
        
        %Q off diagonal elements
        x_sum = 0;
        s0 = 0;
        k = 1;
        for i = 1:n,
            wq = (w(i).*w(i+1:n)) .* Q(k:(k+n-i-1));
            x_sum = x_sum + sum((x(i)-x(i+1:n)).^2 .* wq);
            s0 = s0 + sum(wq);
            k = k + (n-i);
        end
        x_sum = x_sum*2;
        s0 = s0*2;
        
        %Q diagonal elements should be zeros
        pp(j) = 1/2 * (x_sum/s0) / ( sum( (x - ux).^2 .*w )/(sum(w)-1) );

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end

    Vg = 1/( 2*(n+1)*S0.^2 ) * ((2*S1 + S2)*(n-1) - 4*S0.^2);    %gaussian    
    if kk,
        V = Vg;
    else
        fprintf('V<0: %d, V<Vg: %d, V: %d\n',sum(V<0),sum(V<Vg),length(V))
        V(V<Vg) = Vg;
    end
    
    pz = -((pp - EC)./sqrt(V)); %large pp gets small p value
    pz = 2*normcdf(pz);
    p = pz;
    
elseif strcmp(mode, 'moranIt'),

    qpdf = option.pdf;
    k = option.k;
%     if_gaussian = option{4};
    
    %estimate t distribution and probability for each sample
    %t distribution with 1 degree of freedom
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);
    [n,m] = size(X);
    
    %normalize Q
    Q(1:n+1:end) = 0;
    [n,m] = size(X);

    %use hypergeometric test for pairwise distance and expression
    %note that pairwise metrics are not completely independent
    %thresholding expression and distance
    th_Q = quantile(Q(:),1-k/(n-1));
%     idx_Q = Q>th_Q;
%     Q = idx_Q;
    
    EI = -1/(n-1);
    S0 = 2 * sum(Q(:));
    S1 = 2 * sum(Q(:).^2);
    S2 = 4 * sum(sum(Q,1).^2);
%     A = n*((n^2 - 3*n + 3) * S1 - n*S2 + 3*S0^2);
%     B = (n^2 - n) * S1 - 2*n*S2 + 6*S0^2;
%     C = (n-1)*(n-2)*(n-3)*S0^2;
    
    pp = zeros(1,m) + EI;
%     V = ones(1,m);
    
%     Q= squareform(Q);
    
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        
        if sum(x>10) < 1e-10,
            continue
        end        
       
    %     %thresholding method 1
    %     x = double(x>1);
    %     P = pmean_w(x,w,'geomean0');
    %     th_P = 0.5;

        %thresholding method 2
        P = squareform(pmean_w(x,w,'geomean0'));
        P(1:n+1:end) = 0;
        th_P = quantile(x(x>0),0.05);

        P = P>th_P;

%         D = ( sum((x.^4).*w) / sum(w) )/...
%             ( sum((x.^2).*w) / sum(w) ).^2;
%         EI2 = (A - D*B)/C;
%         V(j) = EI2 - EI^2;
        
        x_sum = sum(Q(:) & P(:));
        s0 = sum(Q(:));        
        
        pp(j) = (x_sum/s0) / (sum(x>0) / n);

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    
    %randomization variance sometimes gives erroneous result
    Vg = 1/(S0^2 * (n^2 - 1)) * ( n^2*S1 - n*S2 + 3*S0^2) - EI^2; %gaussian
%     if kk,
%         V = Vg;
%     else
%         fprintf('V<0: %d, V<Vg: %d, V: %d\n',sum(V<0),sum(V<Vg),length(V))
%         V(V<Vg) = Vg;
%     end
    pz = -((pp - EI)./sqrt(Vg));    %large pp gets small p value
    pz = 2*normcdf(pz);
    p = pz;

elseif strcmp(mode,'aple')
    
    qpdf = option.pdf;
    kk = option.prox;  %select for gaussian or randomization
     
    %use APLE
    Q = qpdf( squareform(pdist(ydata,'euclidean')) , 1);

    %normalize Q
    [n,m] = size(X);
    Q(1:n+1:end) = 0;   %not counting self
    ita = 1/sum(Q(:));
    Q = max(Q * ita, realmin);

    pp = zeros(1,m);
    
    Q(1:n+1:end) = 0;
    
    %need to standardized
    Q = bsxfun(@rdivide,Q,sum(Q,2));
    
    lambda = eig(Q);
    figure; plot(lambda,'.');
    Qd = Q'*Q + lambda'*lambda*eye(n)/n;
    Qdd = diag(Qd);
    Qd(1:n+1:end) = 0;
    
    Q = squareform( (Q+Q')/2 );
    Qd = squareform(Qd);
    
    tic
    for j = 1:m,
        x = double(X(:,j)');
        w = double(W(:,j)');
        x = x - mean_w(x',w');
        
        if sum(x>10) < 1e-10,
            continue
        end
                
        %Q off diagonal elements
        x_sum = 0;
        xd_sum = 0;
        s0 = 0;
        sd0 = 0;
        k = 1;
        for i = 1:n,    %jit
            ww = w(i).*w(i+1:n);
            wq = ( ww .* Q(k:(k+n-i-1)) );
            wqd = ( ww .* Qd(k:(k+n-i-1)) );
            xx = (x(i)*x(i+1:n));
            x_sum = x_sum + sum( xx .* wq);
            xd_sum = xd_sum + sum( xx .* wqd);
            s0 = s0 + sum(wq);  %num needs to be normalized
            sd0 = sd0 + sum(ww);
            k = k + (n-i);
        end
        x_sum = x_sum*2;
        xd_sum = xd_sum*2;
        s0 = s0*2;
        sd0 = sd0*2;
        
        %Q diagonal elements
        %note that w.^2 is used for var(x), because it is added to the
        %sum of the off diagonal elements, usually, the var(x)
        %w is used instead of w.^2
        wqd = w.*w .* Qdd';
        xd_sum = xd_sum + sum(x.*x .* wqd);
        sd0 = sd0 + sum(w.*w);
        
        pp(j) = (x_sum/s0) / (xd_sum/sd0);

        if ~rem(j,1000),
            fprintf('j = %d, ',j);
            toc
        end
    end
    p = pp;
    
else
    p = 1;
end

end

%observations are rows
%variables are columns
%calculate range p value based on order statistics on uniform distribution
%ths parameter gp_th needs to be adjusted
function p = ostat_w(dat, gp_th)
[~,m] = size(dat);
p = ones(1,m);
triul = @(t)reshape(t(~diag(ones(1,size(t, 1)))), size(t)-[1 0]);
for j = 1:m,
    %Rousseeuw and Croux method for gaussian estimation
    med = median(dat(:,j));
%     v = sqrt(var(dat(:,j)));
%     v = 1.4826*median(pdist(dat(:,j),'cityblock'));
    v = squareform( pdist(dat(:,j),'cityblock') );
    v = triul(v);
    v = 1.1926*median(median(v));
    
    if ~v,
        continue
    end
%     gp_th = 0.02;
    i = normcdf( -abs(dat(:,j)-med)/v )>gp_th;
    n = sum(i);
    w = max(dat(i,j)) - min(dat(i,j));
    if ~n || ~w
        continue,
    end
    p(j) = w^(n-1)*( n - (n-1)*w );
end
end