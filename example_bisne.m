%{
This example shows how to use bisne to select genes that have
spatial pattern in tSNE embeddings

Yinqing Li
yinqingl@mit.edu
2016

%}

%load necessary functions
import bioma.data.*

fh_s = {
func_util_bisne,
func_util_pdet,
};
for j = 1:length(fh_s),
    fh = fh_s{j};
    fh_str = structvars(fh);
    for i = [1:size(fh_str,1)],
        eval(fh_str(i,:));
    end
end


% %generate some synthetic data
num_genes = 5000;
num_cells_per_group = 50;
num_groups = 4;

i_groups = {};
for i = 1:num_groups,
    i_groups{i} = [(i-1)*num_cells_per_group+1:i*num_cells_per_group];
end
% X = randn(num_genes, length([i_groups{:}]));
% W = ones(size(X));    %probability of detection, for real applications, obtain W from estimate_detection_probability_kde
% 
% %add group specific expression to each group
num_genes_rand = 20;
% for i = 1:num_groups,
%     for j = i:num_groups,
%         gene_i = datasample([1:num_genes], num_genes_rand,'Replace',false);
%         cell_ij = [i_groups{i},i_groups{j}];
%         X(gene_i, cell_ij) = X(gene_i, cell_ij)*randn(1,1) + repmat(randn(length(gene_i),1),1,length(cell_ij));
%     end
% end
% %add noise expression
% num_cells_rand = round(size(X,2)/2);
% for i = 1:10,
%     cell_i = datasample([1:size(X,2)], num_cells_rand,'Replace',false);
%     gene_i = datasample([1:num_genes], num_genes_rand,'Replace',false);
%     X(gene_i, cell_i) = X(gene_i, cell_i)*randn(1,1) + repmat(randn(length(gene_i),1),1,length(cell_i));
% end
% 
% X = log(exp(X)+1);
% 
% X_init = X;
% W_init = W;
% 
% save example_dat X_init W_init

load example_dat

num_cells_rand = round(size(X,2)/2);

X = X_init;
W = W_init;

%perform pca and tsne
[pc_w, ~, pcvars_w] = pca_w(X', W');
figure;semilogy(pcvars_w,'.')
i_pc = [1:3];   %select top PC to be used in tsne embedding
d_wct = pdist_w(X', W', 'cosine', pc_w(:,i_pc)); d_wct1 = corr2d(1-d_wct);
ydata_wct = tsne_matlab(squareform(d_wct1),[],2,[],26,'precomputed'); ydata = ydata_wct;
ydata_pre_bisne = ydata;

%visualize cell group in tsne embedding
label = zeros(length([i_groups{:}]),1);
for i = 1:length(i_groups),
    label(i_groups{i}) = i;
end
f = figure;scatter(ydata(:,1),ydata(:,2),21,label,'filled');

%perform bisne
p_ostat = local_expr_dmap(X',W',ydata,struct('mode','ostat','gp_th',0.02,'th_qtle',1.1));
p_moran = local_expr_dmap(X',W',ydata,struct('mode','moranI','pdf',@(x) tpdf(x,1),'prox',1));

%visualize top scored genes
pp = {prod(p_ostat,2), p_moran};

p = pp{1};
[Y,I] = sort(p);
Xj = X(I(1:50),:);
scatter_series(ydata, Xj)

%search for cut-off score
k = 1; %foreach k in [1,2]
p = pp{k};
[Y,I] = sort(p);    %a ranked list of genes according to bisne statistics
warning('off','stats:statrobustfit:IterationLimit');

ii = [1e1:100]; %test from the first 10 to 50 genes in the ranked list
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
        p = randperm(size(X,2));
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
figure;
for j = 1:pmut_size,
    plot(pcvars_w_diff(1,:) - pcvars_w_diff(1+j,:),'.','Color',get_line_color(j, pmut_size));   %compare to matrix with permutated columns
    hold on
end

pp_I = [55,18]+1e1;   %set cut-off of ranks according to the permutation test
i_glist = [];
for i = 1:length(pp),
    p = pp{i};
    [Y,I] = sort(p);
    i_glist = union(i_glist, I(1:pp_I(i)));
end

X = X_init(i_glist,:);
W = W_init(i_glist,:);

[pc_w, ~, pcvars_w] = pca_w(X', W');
figure;semilogy(pcvars_w,'.')
i_pc = [1:6];   %select top PC to be used in tsne embedding
d_wct = pdist_w(X', W', 'cosine', pc_w(:,i_pc)); d_wct1 = corr2d(1-d_wct);
ydata_wct = tsne_matlab(squareform(d_wct1),[],2,[],26,'precomputed'); ydata = ydata_wct;
ydata_post_bisne = ydata;

%visualize cell group in tsne embedding
label = zeros(length([i_groups{:}]),1);
for i = 1:length(i_groups),
    label(i_groups{i}) = i;
end
f = figure;scatter(ydata(:,1),ydata(:,2),21,label,'filled');

save example_dat_bisne X_init W_init ydata_pre_bisne ydata_post_bisne