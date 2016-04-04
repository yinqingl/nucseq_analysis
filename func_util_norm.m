%{
helper functions for normalization

Yinqing Li
yinqingl@mit.edu
2014

%import functions
fh = func_util_norm;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

Copyright (c) 2015, Yinqing Li
%}
 
function func_h = func_util_norm
func_h = struct();
func_h.norm_DESeq = @norm_DESeq; %normalization using size factor method introduced in DESeq, in house implementation
func_h.norm_rDESeq = @norm_rDESeq; %robust DESeq
func_h.norm_rDESeq_refmat = @norm_rDESeq_refmat; %normalize data to a given reference matrix
func_h.norm_TMM = @norm_TMM; %normalization using tmm, in house implementation
func_h.norm_kde2d = @norm_kde2d; %normalization using kde based method
func_h.norm_kde2d_refmat = @norm_kde2d_refmat; %normalize data to a given reference matrix

%make sure these functions are in the directory
%{
norm_by_maplot_kde_log; %internal helper function for norm_kde2d
lowess_glm; %internal helper function for norm_by_maplot_kde_log
kde2d; %internal helper function for norm_by_maplot_kde_log
%}

end

function [X_norm, S_norm] = norm_DESeq(X)
    k_i = X./(repmat(geomean(X,2),1, size(X,2)));
    k_i_nz = geomean(X,2) ~= 0;
    if sum(k_i_nz) == 0,
        disp('normalization failed')
    end
    for i = [1:size(X,2)],
        S_norm(i) = median(k_i(k_i_nz & (X(:,i)>0),i));
    end
    X_norm = X*diag(1./S_norm);
end

function [X_norm, S_norm] = norm_rDESeq(X_in,if_log)
    if max(max(X_in)) < 100,
        disp('X should be on the linear scale.')
        X_norm = [];
        S_norm = [];
        return
    end
    X_in = double(X_in);
    
    if ~exist('if_log', 'var')
        X = X_in;
        k_i = zeros(size(X));
        for i = [1:size(X,1)],
            k_i(i,X(i,:)~=0) = X(i,X(i,:)~=0)/geomean(X(i,X(i,:)~=0));
        end
        S_norm = zeros(1,size(X,2));    
        for i = [1:size(X,2)],
            S_norm(i) = median(k_i(k_i(:,i)~=0,i));
        end
        X_norm = X*diag(1./S_norm);
    else
        X = log(X_in+1);
        gm = exp(mean(X,2)-1);  %use log(x+1) to approximate geomean
        k_i = zeros(size(X));
        for i = [1:size(X,1)],
            k_i(i,X(i,:)~=0) = X(i,X(i,:)~=0)/gm(X(i,:)~=0);
        end
        S_norm = zeros(1,size(X,2));    
        for i = [1:size(X,2)],
            S_norm(i) = median(k_i(k_i(:,i)~=0,i));
        end        
    end
end

%with reference matrix
function [X_norm, S_norm] = norm_rDESeq_refmat(X, i_X, i_X_ref, n_ref)
    if max(max(X)) < 100,
        disp('X should be on the linear scale.')
        X_norm = [];
        S_norm = [];
        return
    end
    
    j_i = find(sum(X(:,i_X_ref)>0,2)>n_ref);
    disp(sprintf('number of reference genes: %d', length(j_i)));
    k_i = zeros(size(X));
    for ii = 1:length(j_i),
        i = j_i(ii);
        k_i(i,X(i,:)~=0) = X(i,X(i,:)~=0)/geomean(double(X(i,i_X_ref(X(i,i_X_ref)~=0))));
    end
    S_norm = zeros(1,size(X,2));    
    for i = [1:size(X,2)],
        S_norm(i) = median(k_i(j_i(k_i(j_i,i)~=0),i));
    end
    X_norm = X*diag(1./S_norm);
    X_norm = X_norm(:,i_X);
    S_norm = S_norm(i_X);

end

function [X_norm, S_norm] = norm_TMM(X, ref)
    logratioTrim = 0.3;
    sumTrim = 0.05;
    Acutoff = -1e10;
    
    nR = sum(ref);
    ref_nR = ref/nR;
    nR_ref = 1./ref - 1/nR;
    
    S_norm = zeros(1, size(X,2));
    for i = [1:size(X,2)],
        obs = X(:,i);
        nO = sum(obs);
        obs_nO = obs/nO;
        logR = log2(obs_nO) - log2(ref_nR);
        absE = (log2(obs_nO)+log2(ref_nR))/2;
%         logR = log2(obs_nO+1) - log2(ref_nR+1);
%         absE = (log2(obs_nO+1)+log2(ref_nR+1))/2;
        v = (1./obs - 1/nO) + nR_ref;
%         fin = (obs~=0) & (ref~=0) & (absE > Acutoff);
        fin = isfinite(logR) & isfinite(absE) & (absE > Acutoff);
        
        logR = logR(fin);
        absE = absE(fin);
        v = v(fin);
        
        n = sum(fin);
        [loL, hiL] = list(quantile(logR,[logratioTrim, 1-logratioTrim]));
        [loS, hiS] = list(quantile(absE, [sumTrim, 1-sumTrim]));
        keep = (logR>=loL & logR<=hiL) & (absE>=loS & absE<=hiS);
        
        S_norm(i) = 2^(sum(logR(keep)./v(keep)) / sum(1./v(keep)));
%         S_norm(i) = 2^(mean(logR(keep)));
    end
    
    X_norm = X*diag(1./S_norm);
end

function [X_norm, S_norm] = norm_kde2d(X,n_ref)
    if max(max(X)) < 100,
        disp('X should be on the linear scale.')
        X_norm = [];
        S_norm = [];
        return
    end
    
    X_log = log(X+1);
    if ~exist('n_ref', 'var')
        n_ref = 20;
    end
    ma_norm_list = {};
    [dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(X_log);
    I_genes = [1:num_genes];
    [n_genes_detected,idx_sample] = sort(find_num_detected_genes(dat));

    ma_norm_idx = 1;
    if 1/((0.85 - 0.8)/n_ref) > size(X,2),
        qs = linspace(0.8,0.95,n_ref);
    else
        qs = linspace(0.8,0.85,n_ref);
    end

    for i = [1:length(qs)],
        ma_norm_list{ma_norm_idx} = struct();
        q = qs(i)
        ma_norm_list{ma_norm_idx}.quantile = q;
        s_cell_i = idx_sample(find(n_genes_detected>=quantile(n_genes_detected,q),1));
        s_cell_name = cellnames{s_cell_i}
        s_cell = dat(:,s_cell_i); %select the standard cell based on quantile
        [dat_norm_scaling, scaling_factors] = norm_by_maplot_kde_log(dat, I_genes, s_cell, struct('method','scaling'));
        ma_norm_list{ma_norm_idx}.s_cell_i = s_cell_i;
        % ma_norm_list{ma_norm_idx}.dat_norm_scaling = dat_norm_scaling;
        ma_norm_list{ma_norm_idx}.scaling_factors = scaling_factors;
        ma_norm_idx = ma_norm_idx+1;
    end
%     save dat_norm_scaling_20150105_2 ma_norm_list -v7.3

    scaling_factors = ones(length(ma_norm_list), size(ma_norm_list{1}.scaling_factors,2));
    for i = [1:length(ma_norm_list)],
        scaling_factors(i,:) = ma_norm_list{i}.scaling_factors;
    end

%     j = [1,4,6,7,18];
%     scaling_factors(j,:) = [];

    %use Least trimmed squares as regression method
    scaling_factors = log(scaling_factors);
    scaling_factors = scaling_factors - repmat(scaling_factors(:,ma_norm_list{1}.s_cell_i),1,size(scaling_factors,2));
    std_sf = std(scaling_factors,[],1);
    i = std_sf<quantile(std_sf,0.95);

    scaling_factors_a = scaling_factors(:,i);
    a0 = rand(size(scaling_factors_a,1),1)*1e-1;
    a0(1) = 0;
    y = @(a) sum(var(scaling_factors_a - repmat(a,1,size(scaling_factors_a,2)),0,1));
    a = fminsearch(y, a0,struct('MaxFunEvals',5e5,'MaxIter',5e5));
    scaling_factors = scaling_factors - repmat(a,1,size(scaling_factors,2));

    S_norm = exp(median(scaling_factors,1));
    X_norm = (exp(X_log)-1)./repmat(S_norm,size(X_log,1),1);
end

function [X_norm, S_norm] = norm_kde2d_refmat(X, i_X, i_X_ref, n_ref)
    if max(max(X)) < 100,
        disp('X should be on the linear scale.')
        X_norm = [];
        S_norm = [];
        return
    end

    X_log = log(X+1);
    if ~exist('n_ref', 'var')
        n_ref = 20;
    end

    [dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(X_log);
    I_genes = [1:num_genes];
    [n_genes_detected,idx_sample] = sort(find_num_detected_genes( dat(:,i_X_ref) ));

    ma_norm_idx = 1;
    if 1/((0.85 - 0.8)/n_ref) > size( X(:,i_X_ref) ,2)
        n_ref = floor((0.85 - 0.8)*size( X(:,i_X_ref) ,2))
    end
    if 1/((0.85 - 0.8)/n_ref) > size( X(:,i_X_ref) ,2),
        qs = linspace(0.8,0.95,n_ref);
    else
        qs = linspace(0.8,0.85,n_ref);
    end
    
    ma_norm_list = cell(1,length(qs));
    for i = [1:length(qs)],
        ma_norm_list{ma_norm_idx} = struct();
        q = qs(i)
        ma_norm_list{ma_norm_idx}.quantile = q;
        s_cell_i = i_X_ref(idx_sample(find(n_genes_detected>=quantile(n_genes_detected,q),1)));
        s_cell_name = cellnames{s_cell_i}
        s_cell = dat(:,s_cell_i); %select the standard cell based on quantile
        [~, scaling_factors] = norm_by_maplot_kde_log( dat(:,i_X) , I_genes, s_cell, struct('method','scaling'));
        ma_norm_list{ma_norm_idx}.s_cell_i = s_cell_i;
        % ma_norm_list{ma_norm_idx}.dat_norm_scaling = dat_norm_scaling;
        ma_norm_list{ma_norm_idx}.scaling_factors = scaling_factors;
        ma_norm_idx = ma_norm_idx+1;
    end
%     save dat_norm_scaling_20150105_2 ma_norm_list -v7.3

    scaling_factors = ones(length(ma_norm_list), size(ma_norm_list{1}.scaling_factors,2));
    for i = [1:length(ma_norm_list)],
        scaling_factors(i,:) = ma_norm_list{i}.scaling_factors;
    end

%     j = [1,4,6,7,18];
%     scaling_factors(j,:) = [];quit

    %use Least trimmed squares as regression method
    scaling_factors = log(scaling_factors);
    scaling_factors = scaling_factors - repmat(scaling_factors(:,1),1,size(scaling_factors,2));
    std_sf = std(scaling_factors,[],1);
    i = std_sf<quantile(std_sf,0.95);

    scaling_factors_a = scaling_factors(:,i);
    a0 = rand(size(scaling_factors_a,1),1)*1e-1;
    a0(1) = 0;
    y = @(a) sum(var(scaling_factors_a - repmat(a,1,size(scaling_factors_a,2)),0,1));
    a = fminsearch(y, a0,struct('MaxFunEvals',5e5,'MaxIter',5e5));
    scaling_factors = scaling_factors - repmat(a,1,size(scaling_factors,2));

    S_norm = exp(median(scaling_factors,1));
    X_log = X_log(:,i_X);
    X_norm = (exp(X_log)-1)./repmat(S_norm,size(X_log,1),1);
end
