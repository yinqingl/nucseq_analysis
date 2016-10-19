%{
normalize all data points to a reference data point

Yinqing Li
yinqingl@csail.mit.edu
2015

calculate normalization on selected subset of genes
input:
dat, row: genes, col: samples, log(tpm+1)
I_genes, index for selected subset of genes
could be genes with low cv or certain expression level
s_cell, a standard cell to which data is normalized
output:
dat_norm, normalized tpm matrix log(tpm+1)

Copyright (c) 2015, Yinqing Li
%}
function [dat_norm, s] = norm_by_maplot_kde_log(dat, I_genes, s_cell, option)

[n_row,n_col] = size(dat);
dat_norm = zeros(size(dat));
s = ones(1,n_col);

if nargin<3,
    %get a standard cell
    s_cell = median(dat,2);
end

if nargin < 4,
    option = struct();
    option.method = 'lowess';
end

dat_x = s_cell;
for i=[1:n_col],
    disp(sprintf('normalize sample %d',i));
    dat_y = dat(:,i);
    if sum(dat_x - dat_y),
        [dat_norm(:,i),s(i)] = get_dat_norm_on_subset_genes(dat_x, dat_y, I_genes, option);
        %{
        plotname = sprintf('norm_by_maplot_kde_log option = %s %d',option.method,i);
        title(plotname,'Interpreter','none')
        xlabel('a')
        ylabel('m')
        print(gcf,'-r250','-dpng',['analysis_figures' filesep plotname '.png'])
        close(gcf)
        %}
    else
        dat_norm(:,i) = dat_y;
    end
end
end

%return monotonic vector
function [v,m] = get_mnt_vector(v, m)
[v, v_i] = sort(v);
m = m(v_i);
v_diff_i = find(diff(v)==0)+1;
v(v_diff_i) = [];
m(v_diff_i,:) = [];
end

%normalize dat_y to dat_x
%calculate normalization factor based on subset of genes by I_genes
function [dat_y_norm, s] = get_dat_norm_on_subset_genes(dat_x, dat_y, I_genes, option)
log_dat_x = dat_x;
log_dat_y = dat_y;
m = log_dat_y - log_dat_x;
a = (log_dat_y + log_dat_x)/2;
%only choose expressed genes from I_genes
% I_exp = intersect(find(dat_x(I_genes)),find(dat_y(I_genes)));
I_exp = intersect(find(dat_x(I_genes)>1),find(dat_y(I_genes)>1));
% I_exp = intersect(I_exp,find((a_subset > 5 & a_subset < 25)));
m_subset = m(I_genes);m_subset = m_subset(I_exp);
a_subset = a(I_genes);a_subset = a_subset(I_exp);

%use kde to select I_genes
[bandwidth,density,X,Y]=kde2d([a_subset(:),m_subset(:)]);
% figure()
% surf(X,Y,density,'LineStyle','none'), view([0,70])
% colormap hot, hold on, alpha(.8)
% set(gca, 'color', 'blue');
% plot(X(:),Y(:),'b.','MarkerSize',5)

% %only count maxmimum with M within (-2,+2) range
% i = (Y(:)>-3 & Y(:)<3);
i = isfinite(Y(:));
[density_max,density_max_I] = max(density(i));                                    
j = density(:)>=density_max;
density((~i)&j) = density_max*0.95;
i = find(i);
density_max_I = i(density_max_I);
density(density==density_max) = density_max*0.95;
density(density_max_I) = density_max;
density = density/max(density(:));

vq = griddata(X(:),Y(:),density(:),a_subset,m_subset);
% figure()
% plot3(a_subset,m_subset,vq,'.')
level = 0.7;
% vq_th = quantile(vq, level);                                                  %arbitrarily picked 70% quantile
vq_th = level;
% I_density = (vq>vq_th);
% vq = vq(I_density);
% a_subset = a_subset(I_density);
% m_subset = m_subset(I_density);

% remove separated area other than area containing maxmimum peak
% figure()
% imagesc(density)
density_bw = im2bw(density,level);
% figure()
% imagesc(density_bw);
density_CC = bwconncomp(density_bw,4);

if length(density_CC.PixelIdxList)>1,
    for i = [1:length(density_CC.PixelIdxList)],
        %remove areas that do not contain maxmimum peak
        if sum(density_CC.PixelIdxList{i} == density_max_I) == 0,
            density(density_CC.PixelIdxList{i}) = vq_th;
        end
    end
end
density = density - vq_th;
density(density<=0) = 0;
density = density/max(density(:));

% remove local maxmimum peaks
% the local maxmimum peak is defined as
% the distance between the peak point and the point that has higher density
% is larger than a defined radius
while(1)
density_i = find(density>0);
[density_Y, density_I] = sort(density(density_i),'descend');
density_points = [X(:),Y(:)];
density_points = density_points(density_i,:);
density_D = squareform(pdist(density_points,'euclidean'));
density_D_Y = zeros(size(density_I));
for i = [2:length(density_I)],
    density_D_Y(i) = min(density_D(density_I(i),density_I(1:i-1)));
end
% figure;plot((density_D_Y),'.')
density_localmax_I = find(density_D_Y>0.75,1);                             %define the radius of the second peak on average
if ~isempty(density_localmax_I),
    %line search for threshold that separate the two area
    density_localmax_Y = density_Y(density_localmax_I);
    density_localmax_I = density_i(density_I(density_localmax_I));
%     density_ths = fliplr(linspace(0,density_localmax_Y,100));
    density_ths = fliplr(linspace(0,ceil(density_localmax_Y*100),ceil(density_localmax_Y*100)))/100;
    for i = [1:length(density_ths)],
        density_bw = im2bw(density,density_ths(i));
        density_CC = bwconncomp(density_bw,4);
        if length(density_CC.PixelIdxList)>1,
            break
        end
    end
    for i = [i:length(density_ths)],
        density_bw = im2bw(density,density_ths(i));
        density_CC = bwconncomp(density_bw,4);
        if length(density_CC.PixelIdxList)==1,
            break
        end
    end
    density_th = density_ths(i-1);
    density_bw = im2bw(density,density_th);
    density_CC = bwconncomp(density_bw,4);
    for i = [1:length(density_CC.PixelIdxList)],
        %remove areas that do not contain maxmimum peak
        if sum(density_CC.PixelIdxList{i} == density_max_I) == 0,
            density(density_CC.PixelIdxList{i}) = density_th;
        end
    end  
    density = density - density_th;
    density(density<=0) = 0;
    density = density/max(density(:));
else
    break;
end
end

density = density/sum(density(:));
vq = griddata(X(:),Y(:),density(:),a_subset,m_subset);
I_density = (vq>0);
vq = vq(I_density);
a_subset = a_subset(I_density);
m_subset = m_subset(I_density);

if strcmp(option.method,'lowess'),
%select a normalization methods
%lowess
m_subset_fit = smooth(a_subset,m_subset,0.5,'lowess');    %specify method, span = 0.5 is more robust than 0.3
% m_subset_fit = smooth(a_subset,m_subset,0.5,'loess');    %specify method
% m_subset_fit = smooth(a_subset,m_subset,0.5,'rlowess');    %specify method
% m_subset_fit = smooth(a_subset,m_subset,0.5,'rloess');    %specify method
end
if strcmp(option.method,'linear_regression'),
%linear regression
% X = ones(size(m_subset));
% b_fit = robustfit(X,m_subset);
% m_subset_fit = [ones(size(m_subset)),X]*b_fit;
X = [a_subset];
b_fit = robustfit(X,m_subset);
m_subset_fit = [ones(size(m_subset)),X]*b_fit;
end
if strcmp(option.method,'scaling'),
%scaling based on kde
m_subset_fit = repmat(vq(:)'*m_subset(:)/sum(vq(:)),size(m_subset,1),size(m_subset,2));
end

[a_subset_v,m_subset_fit_v] = get_mnt_vector(a_subset,m_subset_fit);
m_fit = interp1(a_subset_v,m_subset_fit_v,a,'nearest','extrap');
s = exp(mean(m_fit));

%normalization in linear space
log_dat_y_fit = log((exp(log_dat_y)-1)./exp(m_fit)+1);
dat_y_norm = log_dat_y_fit;
dat_y_norm(dat_y==0) = 0;
%{
log_dat_y_fit = log_dat_y - m_fit;
dat_y_norm = log_dat_y_fit;
dat_y_norm(dat_y==0) = 0;                                                   %do not normalized 0 count
%}

%{
density(density<=0)=nan;
figure()
set(gcf,'visible','off');
plot(a,m,'.','MarkerSize',1);
hold on
pcolor(X,Y,density)
set(gca, 'color', [1 1 1]);
shading flat
plot(a,m_fit,'-r');
xlabel('A')
ylabel('M')
%}

% figure();
% set(gcf,'visible','off');
% plot(a_subset,m_subset,'.b');
% hold on;
% plot(a_subset,m_subset_fit,'or');
% 
% figure();
% % set(gcf,'visible','off');
% plot(a,m,'.b');
% hold on;
% plot(a,m_fit,'or');

% h = figure();
% plot(a,m,'.')
% hold on
% plot(a_subset,m_subset-m_subset_fit,'r.')
% 
% m = log_dat_y_fit - log_dat_x;
% plot(a, m, 'bo')

end

function [dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(dat_matrix)
dat = double(dat_matrix);
cellnames = get(dat_matrix, 'ColNames');
genes = get(dat_matrix, 'RowNames');
num_samples = get(dat_matrix, 'NCols');
num_genes = get(dat_matrix, 'NRows');
end

function idx_cellstr = find_idx_of_list2_in_list1(cell_list1,cell_list2)
idx_cellstr = {};
for i = [1:length(cell_list2)],
    idx = find(strncmp(cell_list2{i},cell_list1,length(cell_list2{i})));
    idx_cellstr{i} = idx;
end
end