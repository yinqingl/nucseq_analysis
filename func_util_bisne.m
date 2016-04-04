%{
helper functions for dimensionality reduction using PCA, tSNE, and biSNE

Yinqing Li
yinqingl@mit.edu
MIT, 2015

%import functions
import bioma.data.*
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

Copyright (c) 2015, Yinqing Li
%}
 
function func_h = func_util_bisne
func_h = struct();
func_h.corr2d = @corr2d;
% func_h.hogben_corr_m1 = @hogben_corr_m1; @internal helper for corr2d
func_h.eigenvalue_fit = @eigenvalue_fit;
func_h.pca_tsne = @pca_tsne; %example for PCA and tSNE embedding
func_h.bisne_selection = @bisne_selection; %example for selection of bisne threshold

%make sure these functions are in the directory
%{
tsne_matlab; tsne for matlab
d2p_matlab; internal function for tsne_matlab
tsne_p_matlab; internal function for tsne_matlab
local_expr_dmap; calculate significance of genes based on their expression pattern in tsne embedding
func_util_pdet; helper functions supporting PCA with weights
%}

end

%{
example for PCA and tSNE embedding

input:
X: expression data, row: genes, column: samples
W: probability of detection, obtained from estimate_detection_probability_kde

output:
ydata: tsne 2d embedding

procedure:
(1) select top PC to be used in tsne embedding
(2) compute tsne 2d embedding
%}
function pca_tsne
[pc_w, ~, pcvars_w] = pca_w(X', W');
figure;semilogy(pcvars_w,'.')
i_pc = [1:4];   %select top PC to be used in tsne embedding
d_wct = pdist_w(X', W', 'cosine', pc_w(:,i_pc)); d_wct1 = corr2d(1-d_wct);
ydata_wct = tsne_matlab(squareform(d_wct1),[],2,[],26,'precomputed'); ydata = ydata_wct;
end

%{
example for bisne

input:
X: expression data, row: genes, column: samples
W: probability of detection, obtained from estimate_detection_probability_kde
ydata: tsne 2d embedding

output:
i_glist: index of selected genes

procedure:
(1) test significance of localized gene expression in the tsne 2d embedding
(2) rank genes according to their test significance
(3) evaluate cut-off of ranks for each test statistics
%}
function bisne
%*********************************************
%scoring expression pattern in tSNE embedding
p_ostat = local_expr_dmap(X',W',ydata,struct('mode','ostat','gp_th',0.02,'th_qtle',1.1));
p_ostat_h = local_expr_dmap(X',W',tiedrank(ydata),struct('mode','ostat','gp_th',0.02,'th_qtle',1.1));
p_hg = local_expr_dmap(X',W',ydata,struct('mode','hg','pdf',@tpdf,'k',5));
p_moran = local_expr_dmap(X',W',ydata,struct('mode','moranI','pdf',@(x) tpdf(x,1),'prox',1));
p_moran_h = local_expr_dmap(X'>1.1,W',struct('mode','moranI','pdf',@(x) tpdf(x,1),'prox',1));

%*********************************************
%find a cut-off score
pp = {prod(p_ostat_h,2), p_hg, p_moran, p_moran_th_normpdf};

p = pp{k};  %set k to select one of the bisne statistics
p(sum((X>1.1),2)<10) = inf;    %filter on very lowly expressed gene
p(sum((X>1.1),2)>0.8*size(X,2)) = inf;    %filter on widely expressed gene
[Y,I] = sort(p);    %a ranked list of genes according to bisne statistics
warning('off','stats:statrobustfit:IterationLimit');

ii = [1e2:1e3]; %test from the first 100 to 1000 genes in the ranked list
pmut_size = 5;  %number of permutations

pcvars_w_diff = zeros(1+pmut_size,length(ii));
i_glist = I(1:ii(end));
X_t = X(i_glist,:);
W_t = W(i_glist,:);
tic
for i = 1:length(ii),
    [~,~,pcvars_w] = pca_w(X_t(1:ii(i),:)',W_t(1:ii(i),:)');
    pcvars_w = pcvars_w / sum(pcvars_w);
    y = pcvars_w(:);
    [y_fit,y_fit_l] = eigenvalue_fit(y);
    pcvars_w_diff(1,i) =  sum(y(y_fit_l) - y_fit(y_fit_l));
    
    %permutation test
    for j = 1:pmut_size,
        Xp = X_t;
        Wp = W_t;
        p = randperm(length(i_cells));
        Xp(ii(i),:) = X_t(ii(i),p);
        Wp(ii(i),:) = W_t(ii(i),p);
        [~,~,pcvars_w] = pca_w(Xp(1:ii(i),:)',Wp(1:ii(i),:)');
        pcvars_w = pcvars_w / sum(pcvars_w);
        y = pcvars_w(:);
        [y_fit,y_fit_l] = eigenvalue_fit(y);    %find the top variance
        pcvars_w_diff(1+j,i) =  sum(y(y_fit_l) - y_fit(y_fit_l));   %amount of change in the top variance
    end
    
    if ~rem(i,100),
        fprintf('iteration = %d, ',i);
        toc
    end
end
figure; plot(pcvars_w_diff(1,:),'.')
figure;
for j = 1:pmut_size,
    plot(pcvars_w_diff(1,:) - pcvars_w_diff(1+j,:),'.','Color',get_line_color(j, pmut_size));   %compare to matrix with permutated columns
    hold on
end
figure;
plot(smooth(prod(bsxfun(@minus, pcvars_w_diff(1,:), pcvars_w_diff(2:end,:)) > 0, 1),10),'.')

%*********************************************
%take union of genes
pp = {p_ostat, p_hg, p_moran};
pp_I = [1,1.6e2+132,1.6e2+198,1.6e2+262];   %set cut-off of ranks according to the permutation test
% i_glist = [1:length(genes)];
i_glist = [];
X = dat_matrix_log_norm2(:,i_cells);
W = p_det(:,i_cells);
for i = 1:length(pp),
    p = pp{i};
    p(sum((X>1.1),2)<10) = inf;    %filter on very lowly expressed gene
    p(sum((X>1.1),2)>0.8*size(X,2)) = inf;    %filter on widely expressed gene
    [Y,I] = sort(p);
    i_glist = union(i_glist, I(1:pp_I(i)));
end

end

%{
convert correlation to euclidean distance
assume gaussian noise model

input:
rho: correlation
Var_X: variance of gene expression, for example, Var_X = var(X,[],2);
n: number of samples, for example n = size(X,2);

output:
d: euclidean distance
%}
function d = corr2d(rho, Var_X, n)
%large n approximation for cosine
d = min( sqrt(max(1./max(rho,realmin).^2,1)-1), 1/sqrt(realmin) );

%use hogben method for correlation
tic
if exist('Var_X', 'var')
    if min(size(rho)) == 1,
        disp('rho needs to be a square matrix.\n');
        return
    end
    m = size(rho,1);
    d = ones(size(rho));
    sum_X2 = sqrt(Var_X*(n-1));
    for i = 1:m,
        for j = 1:m,
            sumX2 = sum_X2(i);
            y = @(x) (hogben_corr_m1(x,sumX2,n-2) - max(rho(i,j),1e-6)^2).^2;
            d(i,j) = fminsearch(y, 1);
        end
        if rem(i,10) == 0,
            disp(sprintf('iteration %d',i));
            toc
        end        
    end
end
end

%{
use hogben method to estimate e
%}
function m = hogben_corr_m1(sigma, sumX2, n)
theta = sumX2/(sigma+realmin);
f1 = @(u) (1-(u/theta).^2).^(n/2).*exp(-u.^2/2);
m = (2/pi)^(1/2)*integral(f1,0,theta);
end

%{
evaluate fitting of higher order eigen values using approximated Marcenko-Pastur distribution
calculate significance of deviations of eigenvalues from the fitting

input:
y: eigenvalues

output:
y_fit: fitted eigenvalues
y_fit_l: index of eigenvalues that deviate from fitting by 0.01 significance level
%}
function [y_fit,varargout] = eigenvalue_fit(y)
y = y/sum(y);
x = [1:length(y)]'/length(y);
xp = [ones(size(x)),x.^(1/2)];
warning('off','stats:statrobustfit:IterationLimit');
b = robustfit(xp,y,[],[],'off');
y_fit = xp*b;
z = y-y_fit;
v = squareform( pdist(z,'cityblock') );
triul = @(t)reshape(t(~diag(ones(1,size(t, 1)))), size(t)-[1 0]);
v = triul(v);
v = 1.1926*median(median(v));
z = (z - median(z))/v;
pz = normcdf(-abs(z));
y_fit_l = [1:length(pz) - find(flipud(pz(:)) > 0.01,1) + 1];
pz_i = (pz > 0.01);
b = robustfit(xp(pz_i,:),y(pz_i),[],[],'off');
y_fit = xp*b;
nout = max(nargout,1) - 1;
for k = 1:nout
   varargout{k} = y_fit_l;
end
end
