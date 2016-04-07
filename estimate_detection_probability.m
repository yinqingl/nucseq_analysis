%estimate the probability of detecting a true tpm
%Yinqing Li
%yinqingl@mit.edu
%input: dat, row: genes, col: samples, tpm
%output: p_det, other probability estimation
%to get strawman lambda_k_est
%lambda_k_est = mean(prob_est.lambda_ik_est,2)
%p_det = [1:n_genes]
%p_det = lambda_k_est(prob_est.b_j)


function [p_det,varargout] = estimate_detection_probability(dat,option)
t = [];
tic;
%dat format, rows: genes, columns: samples
%dat should not contain zero row

is_median_l_e_ij = false;
is_bin_k_and_below_lambda_ij = false;
is_show_iteration = false;
is_quantile = true;
is_consolidate_bins = true;
is_p_e_ij_uniform_prior = true;
n_bin = 21;
norm_tol_th_iter = 1e-2;
n_delta_diff_p_e_ij = 5;
max_iter = 200;

if nargin >1,
    option_str = structvars(option);
    for i = [1:size(option_str,1)],
        eval(option_str(i,:));
    end
end

% load tmp_dat
[nrow, ncol] = size(dat);

%estimate initial p_e_ij
d_ij = dat>0;
if is_p_e_ij_uniform_prior,
    epsilon = 1e-2;
    p_e_ij_1 = ones(size(d_ij));
    p_e_ij_1(d_ij==0) = 1-epsilon;
    p_e_ij_0 = 1-p_e_ij_1;
else
    p_e_j_1 = sum(d_ij,2)./size(d_ij,2);   %prob e_j = 1

    p_e_ij_1 = zeros(size(dat));  %prob e_ij = 1
    p_e_ij_1(d_ij==1) = 1;
    for j = [1:nrow],
        p_e_ij_1(j,d_ij(j,:)==0) = p_e_j_1(j);
    end
    p_e_ij_0 = 1-p_e_ij_1;
end

%estimate initial l_e_ij
%if l_e_ij has large variance, then consider more robust estimate
if ~is_median_l_e_ij,
    l_e_ij = dat;
    l_e_j = sum(l_e_ij.*p_e_ij_1,2)./sum(p_e_ij_1,2);
    l_e_j(isnan(l_e_j)) = 0;
    l_e_ij = l_e_j;
else,
%use median to estimate l_e_ij
    l_e_j = zeros(nrow,1);
    for j = [1:nrow],
        l_e_j(j) = median(dat(j,d_ij(j,:)));
    end
    l_e_j(isnan(l_e_j)) = 0;
    l_e_ij = l_e_j;
end

%divide l_e_ij to bins
if is_consolidate_bins,
    n_bin_init = n_bin;
end
if ~is_quantile,
    epsilon = 1;
    [n,bin_log_l_e_j] = hist(log(l_e_j+epsilon),n_bin);
    bin_size = bin_log_l_e_j(2)-bin_log_l_e_j(1);
    bin_l_e_j = [-inf, exp(bin_log_l_e_j(1:end-1)+bin_size/2)-epsilon, inf];
    [n,b_j]= histc(l_e_j,bin_l_e_j);
else,
    [y,l_e_j_i] = sort(l_e_j);
    
    num_nonzero_l_e_j = sum(l_e_j>0);
    num_zero_l_e_j = length(l_e_j) - num_nonzero_l_e_j;
    if num_zero_l_e_j ~= 0,
        n_bin_nonzero = n_bin-1;
    else,
        n_bin_nonzero = n_bin;
    end
    
    b_j = zeros(num_nonzero_l_e_j,1);
    quantile_edge = round(quantile([1:num_nonzero_l_e_j],linspace(0,1,n_bin_nonzero+1)));
    quantile_edge = quantile_edge(:);
    quantile_edges = [quantile_edge(1:end-1),quantile_edge(2:end)];
    for k = [1:size(quantile_edges,1)],
        b_j(quantile_edges(k,1):quantile_edges(k,2)) = k+1;
    end
    if num_zero_l_e_j ~=0,
        b_j = [ones(num_zero_l_e_j,1);b_j];
    else,
        b_j = b_j - 1;
    end
    b_j(l_e_j_i) = b_j;
    
    if is_consolidate_bins,
        epsilon = 1;
        [n,bin_log_l_e_j] = hist(log(l_e_j+epsilon),n_bin);
        bin_size = bin_log_l_e_j(2)-bin_log_l_e_j(1);
        bin_l_e_j = [-inf, exp(bin_log_l_e_j(1:end-1)+bin_size/2)-epsilon, inf];
        [n,b_j_histc]= histc(l_e_j,bin_l_e_j);
        [b_j, n_bin] = consolidate_bins(n_bin, b_j, b_j_histc);
    end
    
    log_l_e_j = log(l_e_j+epsilon);
    bin_log_l_e_j = zeros(n_bin,1);
    for k = [1:n_bin],
        bin_log_l_e_j(k) = mean(log_l_e_j(b_j == k));
    end
end

%estimate initial p_d_ij_1_e_ij_1, e_ij = 1
lambda_ik = ones(n_bin,ncol);
lambda_ik_est = ones(n_bin,ncol);
for k = [1:n_bin],
    if ~is_bin_k_and_below_lambda_ij,
        epsilon = 1e-3;
        lambda_ik(k,:) = sum(d_ij(b_j == k,:),1)/size(d_ij(b_j == k,:),1);
        lambda_ik(k,(lambda_ik(k,:)<epsilon)) = epsilon;
        lambda_ik_est(k,:) = lambda_ik(k,:);
        
        for i = [1:ncol],
            f = @(x) -sum(d_ij(b_j == k,i))*log(x)-sum(log(1-p_e_ij_1(b_j == k,i)*x).*(d_ij(b_j == k,i)==0));
            g = @(x) -sum(d_ij(b_j == k,i))/x+sum(p_e_ij_1(b_j == k,i)./(1-p_e_ij_1(b_j == k,i)*x).*(d_ij(b_j == k,i)==0));
            f_g = @(x) deal(f(x),g(x));
%             %use fzero
%             lambda_ik(k,i) = fzero(g,lambda_ikk(k,i));
            
%             %use fminsearch
%             lambda_ik(k,i) = fminsearch(f,lambda_ik_est(k,i),optimset('MaxFunEvals',1e3,'Display','off'));
            
            %use fmincon
            epsilon = 1e-3;
            lambda_ik(k,i) = fmincon(f_g,lambda_ik_est(k,i),[],[],[],[],epsilon,1-epsilon,[],optimoptions('fmincon','GradObj','on','Display','notify','UseParallel','always'));
        end
        
    else,
        d_ij_k = [];
        for kk = [1:k],
            d_ij_k = [d_ij_k;d_ij(b_j == kk,:)];
        end
        lambda_ik(k,:) = sum(d_ij_k,1)/size(d_ij_k,1);
    end
end
epsilon = 1e-3;
lambda_ik(lambda_ik<=0) = epsilon;
lambda_ik(lambda_ik>=1) = 1-epsilon;
       
p_d_ij_1_e_ij_0 = zeros(size(dat));   %prob d_ij = 1 when e_ij = 0
p_d_ij_0_e_ij_0 = 1-p_d_ij_1_e_ij_0;   %prob d_ij = 0 when e_ij = 0
p_d_ij_1_e_ij_1 = lambda_ik(b_j,:);   %prob d_ij = 1 when e_ij = 1
p_d_ij_0_e_ij_1 = 1-p_d_ij_1_e_ij_1;   %prob d_ij = 0 when e_ij = 1

%estimate p_e_ij_d_ij
p_e_ij_1_d_ij = zeros(size(dat));
p_e_ij_1_d_ij = p_d_ij_0_e_ij_1.*p_e_ij_1./(p_d_ij_0_e_ij_0.*p_e_ij_0 + p_d_ij_0_e_ij_1.*p_e_ij_1);
p_e_ij_1_d_ij(d_ij == 1) = 1;

%update p_e_ij
epsilon = 1e-3;
p_e_j_1 = mean(p_e_ij_1_d_ij,2);   %prob e_j = 1
p_e_j_1(p_e_j_1>=1-epsilon) = 1-epsilon;
p_e_j_1(p_e_j_1<=epsilon) = epsilon;
p_e_ij_1 = zeros(size(dat));  %prob e_ij = 1
p_e_ij_1(d_ij==1) = 1;
for j = [1:nrow],
    p_e_ij_1(j,d_ij(j,:)==0) = p_e_j_1(j);
end
p_e_ij_0 = 1-p_e_ij_1;
p_e_ij_1_d_ij_previter = p_e_ij_1_d_ij;

%iterate
if_iter = true;
i_iter = 0;
diff_p_e_ij_v = [];
while(if_iter),
    
    %estimate initial l_e_ij
    %if l_e_ij has large variance, then consider more robust estimate
    if ~is_median_l_e_ij,
        l_e_ij = dat;
        l_e_j = sum(l_e_ij.*p_e_ij_1,2)./sum(p_e_ij_1,2);
        l_e_j(isnan(l_e_j)) = 0;
        l_e_ij = l_e_j;
    else,
        %use median to estimate l_e_ij
        l_e_j = zeros(nrow,1);
        for j = [1:nrow],
            l_e_j(j) = median(dat(j,d_ij(j,:)));
        end
        l_e_j(isnan(l_e_j)) = 0;
        l_e_ij = l_e_j;
    end

    %divide l_e_ij to bins
    if is_consolidate_bins,
        n_bin = n_bin_init;
    end    
    if ~is_quantile,
        epsilon = 1;
        [n,bin_log_l_e_j] = hist(log(l_e_j+epsilon),n_bin);
        bin_size = bin_log_l_e_j(2)-bin_log_l_e_j(1);
        bin_l_e_j = [-inf, exp(bin_log_l_e_j(1:end-1)+bin_size/2)-epsilon, inf];
        [n,b_j]= histc(l_e_j,bin_l_e_j);
    else,
        [y,l_e_j_i] = sort(l_e_j);
        
        num_nonzero_l_e_j = sum(l_e_j>0);
        num_zero_l_e_j = length(l_e_j) - num_nonzero_l_e_j;
        if num_zero_l_e_j ~= 0,
            n_bin_nonzero = n_bin-1;
        else,
            n_bin_nonzero = n_bin;
        end        

        b_j = zeros(num_nonzero_l_e_j,1);
        quantile_edge = round(quantile([1:num_nonzero_l_e_j],linspace(0,1,n_bin_nonzero+1)));
        quantile_edge = quantile_edge(:);
        quantile_edges = [quantile_edge(1:end-1),quantile_edge(2:end)];
        for k = [1:size(quantile_edges,1)],
            b_j(quantile_edges(k,1):quantile_edges(k,2)) = k+1;
        end
        if num_zero_l_e_j ~=0,
            b_j = [ones(num_zero_l_e_j,1);b_j];
        else,
            b_j = b_j - 1;
        end
        b_j(l_e_j_i) = b_j;
    
        if is_consolidate_bins,
            epsilon = 1;
            [n,bin_log_l_e_j] = hist(log(l_e_j+epsilon),n_bin);
            bin_size = bin_log_l_e_j(2)-bin_log_l_e_j(1);
            bin_l_e_j = [-inf, exp(bin_log_l_e_j(1:end-1)+bin_size/2)-epsilon, inf];
            [n,b_j_histc]= histc(l_e_j,bin_l_e_j);
            [b_j, n_bin] = consolidate_bins(n_bin, b_j, b_j_histc);
        end

        log_l_e_j = log(l_e_j+epsilon);
        bin_log_l_e_j = zeros(n_bin,1);        
        for k = [1:n_bin],
            bin_log_l_e_j(k) = mean(log_l_e_j(b_j == k));
        end
    end
    toc
    %estimate lambda_ik
    lambda_ik = ones(n_bin,ncol);
    lambda_ik_est = ones(n_bin,ncol);
    for k = [1:n_bin],
        if ~is_bin_k_and_below_lambda_ij,
            epsilon = 1e-3;
            lambda_ik(k,:) = sum(d_ij(b_j == k,:),1)/size(d_ij(b_j == k,:),1);
            lambda_ik(k,(lambda_ik(k,:)<epsilon)) = epsilon;
            lambda_ik_est(k,:) = lambda_ik(k,:);
            
            for i = [1:ncol],
                f = @(x) -sum(d_ij(b_j == k,i))*log(x)-sum(log(1-p_e_ij_1(b_j == k,i)*x).*(d_ij(b_j == k,i)==0));
                g = @(x) -sum(d_ij(b_j == k,i))/x+sum(p_e_ij_1(b_j == k,i)./(1-p_e_ij_1(b_j == k,i)*x).*(d_ij(b_j == k,i)==0));
                f_g = @(x) deal(f(x),g(x));
    %             %use fzero
    %             lambda_ik(k,i) = fzero(g,lambda_ikk(k,i));

    %             %use fminsearch
    %             lambda_ik(k,i) = fminsearch(f,lambda_ik_est(k,i),optimset('MaxFunEvals',1e3,'Display','off'));

                %use fmincon
                epsilon = 1e-3;
                lambda_ik(k,i) = fmincon(f_g,lambda_ik_est(k,i),[],[],[],[],epsilon,1-epsilon,[],optimoptions('fmincon','GradObj','on','Display','notify','UseParallel','always'));
            end

        else,
            d_ij_k = [];
            for kk = [1:k],
                d_ij_k = [d_ij_k;d_ij(b_j == kk,:)];
            end
            lambda_ik(k,:) = sum(d_ij_k,1)/size(d_ij_k,1);
        end
    end
    epsilon = 1e-3;
    lambda_ik(lambda_ik<=0) = epsilon;
    lambda_ik(lambda_ik>=1) = 1-epsilon;
    toc
    p_d_ij_1_e_ij_0 = zeros(size(dat));   %prob d_ij = 1 when e_ij = 0
    p_d_ij_0_e_ij_0 = 1-p_d_ij_1_e_ij_0;   %prob d_ij = 0 when e_ij = 0
    p_d_ij_1_e_ij_1 = lambda_ik(b_j,:);   %prob d_ij = 1 when e_ij = 1
    p_d_ij_0_e_ij_1 = 1-p_d_ij_1_e_ij_1;   %prob d_ij = 0 when e_ij = 1

    %estimate p_e_ij_d_ij
    p_e_ij_1_d_ij = zeros(size(dat));
    p_e_ij_1_d_ij = p_d_ij_0_e_ij_1.*p_e_ij_1./(p_d_ij_0_e_ij_0.*p_e_ij_0 + p_d_ij_0_e_ij_1.*p_e_ij_1);
    p_e_ij_1_d_ij(d_ij == 1) = 1;

    %iteration condition
    i_iter = i_iter + 1;
    diff_p_e_ij = norm(p_e_ij_1_d_ij_previter - p_e_ij_1_d_ij);
    
    if length(diff_p_e_ij_v) > n_delta_diff_p_e_ij,
        delta_diff_p_e_ij = mean(abs(diff(diff_p_e_ij_v(end-n_delta_diff_p_e_ij+1:end))));
    else
        delta_diff_p_e_ij = 0;
    end
    
    if diff_p_e_ij < norm_tol_th_iter,
        if_iter = false;
    elseif delta_diff_p_e_ij < norm_tol_th_iter && length(diff_p_e_ij_v) > n_delta_diff_p_e_ij,
        if_iter = false
    elseif i_iter>max_iter,
        if_iter = false;
    else,   %continue iteration
        if is_show_iteration,
            disp(sprintf('iteration: %d',i_iter))
            disp(sprintf('diff_p_e_ij = %f',diff_p_e_ij))
            disp(sprintf('delta_diff_p_e_ij = %f',delta_diff_p_e_ij))
        end
        
        %update p_e_ij
        epsilon = 1e-3;
        p_e_j_1 = mean(p_e_ij_1_d_ij,2);   %prob e_j = 1
        p_e_j_1(p_e_j_1>=1-epsilon) = 1-epsilon;
        p_e_j_1(p_e_j_1<=epsilon) = epsilon;

        p_e_ij_1 = zeros(size(dat));  %prob e_ij = 1
        p_e_ij_1(d_ij==1) = 1;
        for j = [1:nrow],
            p_e_ij_1(j,d_ij(j,:)==0) = p_e_j_1(j);
        end
        p_e_ij_0 = 1-p_e_ij_1;
        p_e_ij_1_d_ij_previter = p_e_ij_1_d_ij;
        diff_p_e_ij_v(end+1) = diff_p_e_ij;
    end
end

%p_det is defined as the probability of e_ij == d_ij
p_det = zeros(size(dat));
p_det(d_ij==1) = 1;
p_e_ij_0_d_ij = 1-p_e_ij_1_d_ij;
p_det(d_ij==0) = p_e_ij_0_d_ij(d_ij==0);
epsilon = 1e-4;
p_det((l_e_j == 0),:) = 1-epsilon;

%show iteration
%{
format shortg
c = clock;
c = strrep(num2str(fix(c)),' ','');

figure();
plot(diff_p_e_ij_v,'-.');
plotname = sprintf('lambda_iteration_med=%d_qtl=%d_bin=%d_iter=%d_%s',is_median_l_e_ij,is_quantile,n_bin,max_iter,c);
title(plotname,'Interpreter','none')
xlabel('iteration')
ylabel('diff_p_e_ij')
print(gcf,'-r250','-dpng',['analysis_figures/' plotname '.png'])

%show lambda
figure();
cmap = lines(ncol+5);
for i = [1:ncol],
    plot(bin_log_l_e_j,lambda_ik(:,i),'-.','color',cmap(i,:))
    hold on
end
plotname = sprintf('lambda_med=%d_qtl=%d_bin=%d_iter=%d_%s',is_median_l_e_ij,is_quantile,n_bin,max_iter,c);
title(plotname,'Interpreter','none')
xlabel('bin_log_l_e_j','Interpreter','none')
ylabel('lambda_ik','Interpreter','none')
print(gcf,'-r250','-dpng',['analysis_figures/' plotname '.png'])
%}

prob_est = struct();
prob_est.b_j = b_j;
prob_est.lambda_ik = lambda_ik;
prob_est.lambda_ik_est = lambda_ik_est;
prob_est.p_e_ij_1 = p_e_ij_1;
prob_est.p_e_ij_1_d_ij = p_e_ij_1_d_ij;
varargout{1} = prob_est;

% figure()
% imagesc(log(dat(find(sum(dat,2)),:))+1);
% figure()
% imagesc(p_det(find(sum(dat,2)),:));
end

function [b_j_consolidate, n_bin_unique] = consolidate_bins(n_bin, b_j, b_j_histc)
b_2_b_histc = zeros(1,n_bin);
for k = [1:n_bin],
    j_index = (b_j == k);
    b_j_index = b_j_histc(j_index);
    if length(unique(b_j_index)) == 1,
        b_2_b_histc(k) = unique(b_j_index);
    end
end
b_2_consolidate_b = [1:n_bin];
b_2_consolidate_b(1) = 1;
n_bin_unique = 1;
for k = [2:n_bin],
    if b_2_b_histc(k) ~= 0,
        if b_2_b_histc(k) == b_2_b_histc(k-1),
            b_2_consolidate_b(k) = n_bin_unique;
        else,
            n_bin_unique = n_bin_unique + 1;
            b_2_consolidate_b(k) = n_bin_unique;
        end
    else,
        n_bin_unique = n_bin_unique + 1;
        b_2_consolidate_b(k) = n_bin_unique;
    end
end
b_j_consolidate = size(b_j);
for k = [1:n_bin],
    b_j_consolidate(b_j == k) = b_2_consolidate_b(k);
end
b_j_consolidate = b_j_consolidate(:);
end