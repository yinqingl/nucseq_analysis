%{
helper functions for fitting a path along data points

Yinqing Li
yinqingl@mit.edu
2015

%import functions
fh = func_util_pathfit;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

Copyright (c) 2015, Yinqing Li
%}

function func_h = func_util_pathfit
func_h = struct();
% func_h.ginput_interactive = @ginput_interactive; %internal helper function

func_h.find_st = @find_st; %example for fitting a linear structure using LKH or mst
func_h.find_st_bspline = @find_st_bspline; %example for bspline fitting for tsne 2d embedding

%make sure these functions are in the directory
%{
LKH_tsp; %index cells using LKH algorithm
bspline_2dfit; %fitting bspline to 2d embedding
%}

end

%{
find a linear structure using LKH or mst

input:
ydata: n by 2 tsne embedding matrix
X: m by n, gene expression matrix, m genes, n cells

output:
st: order of cells along the path

procedure
(1) select outliers to be removed using polygon selection tool
(2) select the start and the end of the path using polygon selection tool
(3) use LKH algorithm to find the path
%}
function find_st
%remove outlier
figure;scatter(ydata(:,1),ydata(:,2),21,sum(X>1.1,1),'filled');
[xv,yv] = ginput_interactive();   %select polygon
i = find(inpolygon(ydata(:,1),ydata(:,2),xv,yv));
ydata(i,:) = [];

%organize data in upper triangle matrix
d = pdist(ydata,'euclidean');
d = ceil(d/min(d)*10);
%select beginning and end
figure;scatter(ydata(:,1),ydata(:,2),21,sum(X>1.1,1),'filled');
title('select a cell as the start of the path');
[xv,yv] = ginput_interactive();   %select polygon
i_begin = find(inpolygon(ydata(:,1),ydata(:,2),xv,yv));
title('select a cell as the end of the path');
[xv,yv] = ginput_interactive();   %select polygon
i_end = find(inpolygon(ydata(:,1),ydata(:,2),xv,yv));
%use LKH
st = LKH_tsp(d, i_begin, i_end);

figure;scatter(ydata(:,1),ydata(:,2),21,st,'filled');

%use mst
ug_W = squareform(pdist(ydata,'euclidean'));
k = 6+1;    %knn graph
for i = 1:size(ug_W,1),
    [a,b] = sort(ug_W(i,:));
    ug_W(i,a>k) = 0;
end
ug_W = sparse(ug_W);
[ST,pred] = graphminspantree(ug_W,1);

figure; scatter(ydata(:,1),ydata(:,2),21,sum(X>1.1,1),'filled');
hold on
[i,j,s] = find(ST);
for k = 1:length(i);
    x = [ydata(i(k),1);ydata(j(k),1)];
    y = [ydata(i(k),2);ydata(j(k),2)];
    line(x,y,...
        'LineWidth',1,...
        'Color',[1 0 0]);
end
end

%{
find a linear structure using bspline

input:
ydata: n by 2 tsne embedding matrix
X: m by n, gene expression matrix, m genes, n cells

parameters:
d_th: threshold for epsilon graph
k: order of bspline
option.lambda: ratio between bspline total length and orthogonal variance

output:
t_hat: pseudotime of cells along the path

procedure
(1) select the start and the end of the path using polygon selection tool
(2) use bspline fitting to find the path
%}
function find_st_bspline
%use shortest distance to initialize bspline
figure; scatter(ydata(:,1),ydata(:,2),21,sum(X>1.1,1),'filled');

%find affinity graph
ug_W = squareform(pdist(ydata,'euclidean'));

%sparsify graph
%{
%use knn
k = 6+1;
for i = 1:size(ug_W,1),
    [a,b] = sort(ug_W(i,:));
    ug_W(i,b(k+1:end)) = 0;
end
%}

%use limit distance
d_th = 3;   %epsilon graph
% d_th = quantile(ug_W(:),1 - 6/(size(ydata,1)-1));
for i = 1:size(ug_W,1),
    ug_W(i,ug_W(i,:)>d_th) = 0;
end

ug_W = sparse(ug_W);

%select path start and end points
[xv,yv] = ginput_interactive();   %select polygon
i_begin = find(inpolygon(ydata(:,1),ydata(:,2),xv,yv));
[xv,yv] = ginput_interactive();   %select polygon
i_end = find(inpolygon(ydata(:,1),ydata(:,2),xv,yv));
[dist, path, pred] = graphshortestpath(ug_W, i_begin, i_end);

%smooth path
d_th = 3;
ydata_path = zeros(length(path), 2);
for i = 1:length(path),
    ii = [path(i), find(ug_W(path(i),:) > 0 & ug_W(path(i),:) < d_th)];
    ydata_path(i,:) = mean(ydata(ii,:),1);
end

figure; scatter(ydata(:,1),ydata(:,2),21,label(i_cells),'filled');
hold on
for i = 1:length(path)-1;
    x = [ydata_path(i,1);ydata_path(i+1,1)];
    y = [ydata_path(i,2);ydata_path(i+1,2)];
    line(x,y,...
        'LineWidth',1,...
        'Color',[1 0 0]);
end
scatter(ydata_path(:,1),ydata_path(:,2),21,'filled');

%construct initial spline
n = length(ydata_path);
k = 3+1;
m_knot = augknt(linspace(0,1,n-2),k);
% m_knot = optknt(ydata_path(:,1)',4);
m_cv = ydata_path';
sp = spmak(m_knot,m_cv);
figure;
scatter(ydata(:,1),ydata(:,2),21,label(i_cells),'filled');
hold on
fnplt(sp);
scatter(ydata_path(:,1),ydata_path(:,2),21,'filled');

%fit bspline
option = {};
option.is_show_iteration = true;
option.is_show_f1_f2 = true;
option.norm_tol_th_iter = 1e-2;
option.n_delta_diff_f = 3;
option.max_iter = 300;
option.lambda = 5;
option.mode = 'tlen';
option.fixedse = true;
[sp_fit,t_hat] = bspline_2d_fitting(ydata,sp,option);

figure;
scatter(ydata(:,1),ydata(:,2),21,t_hat,'filled');
hold on
fnplt(sp_fit);
scatter(ydata_path(:,1),ydata_path(:,2),21,'filled');
end

%{
select data points using polygon

example:
[xv,yv] = ginput_interactive();   %select polygon
in = inpolygon(ydata(:,1),ydata(:,2),xv,yv);
%}
function [xs,ys] = ginput_interactive()
xs = [];
ys = [];
while 1,
    [x,y] = ginput(1);
    if length(x),
        xs(end+1) = x;
        ys(end+1) = y;
        hold on;
        plot(xs,ys);
    else,
        break
    end
end

end

%{
use LKH to solve tsp

input:
d: pairwise distance
i_begin: beginning of a path
i_end: end of a path
if i_end is not given, i_end = i_begin

output:
st: indexing of data points along the path

need LKH-2.0.2\LKH in system path
%}
function st = LKH_tsp(d,i_begin,i_end)
if ~exist('i_end', 'var')
    i_end = i_begin;
end

if size(d,1) ~= 1,
    d = squareform(d);
    d = ceil(d/min(d)*10);
    d = squareform(d);
else
    d = squareform(d);
end
n = size(d,1);
% d(i_begin,i_end) = min(min(d(:)),1);
d(i_begin,i_end) = 0;
d(i_end,i_begin) = d(i_begin,i_end);
% %keep only knn
% k = 50;
% if n > k,
%     for i = 1:n,
%         [a,b] = sort(d(i,:));
%         d(i,b(k+1:end)) = a(end);
%     end
% end
d = squareform(d);

addpath('LKH-2.0.2')
rng('shuffle');
% fid = ceil(rand*1e6+1);
fid = 1;
filename_tsp = ['LKH-2.0.2' filesep sprintf('st%d.tsp',fid)];
filename_par = ['LKH-2.0.2' filesep sprintf('st%d.par',fid)];
filename_out = ['LKH-2.0.2' filesep sprintf('st%d_tour.txt',fid)];

%write inputs
fh = fopen(filename_tsp,'w');
fprintf(fh,'NAME: tsp\n');
fprintf(fh,'TYPE: TSP\n');
fprintf(fh,'COMMENT: Matlab\n');
fprintf(fh,'DIMENSION: %d\n',n);
fprintf(fh,'EDGE_WEIGHT_TYPE: EXPLICIT\n');
fprintf(fh,'EDGE_WEIGHT_FORMAT: UPPER_ROW\n');
fprintf(fh,'EDGE_WEIGHT_SECTION\n');
% for i = 1:length(d)
%     fprintf(fh,'%d',d(i));
%     if rem(i,10),
%         fprintf(fh,' ');
%     else
%         fprintf(fh,'\n');
%     end
% end
n = length(d);
a = rem(n,10);
b = (n-a)/10;
c = num2str(reshape(d(1:n-a),10,b)','%d ');
for i = 1:b,
    fprintf(fh, '%s\n', c(i,:));
end
c = num2str(d(n-a+1:n),'%d ');
fprintf(fh, '%s\n', c);
fprintf(fh,'EOF\n');
fclose(fh);

%write LKH parameters
fh = fopen(filename_par,'w');
fprintf(fh,'PROBLEM_FILE= %s\n',filename_tsp);
fprintf(fh,'OPTIMUM 378032\n');
fprintf(fh,'MOVE_TYPE = 5\n');
fprintf(fh,'PATCHING_C = 3\n');
fprintf(fh,'PATCHING_A = 2\n');
fprintf(fh,'RUNS = 10\n');
fprintf(fh,'OUTPUT_TOUR_FILE = %s\n',filename_out);
fclose(fh);

cmd = ['LKH-2.0.2' filesep 'LKH ' filename_par];
tic, system(cmd); toc

%load LKH tour
st = importdata(filename_out);
st = st.textdata;
for i = 1:size(st,1),
    if strmatch('TOUR_SECTION',st{i,1}),
        i_seq_begin = i+1;
    elseif strmatch('-1',st{i,1}),
        i_seq_end = i-1;
    end
end
st = cellfun(@str2num,st([i_seq_begin:i_seq_end],1));

%shift st so that it begins with i_begin
i = find(st==i_begin);
j = find(st==i_end);
if i<j,
    st = flipud(st);
end
[~,st] = sort(st);
st = mod(st + (length(st) - st(i_begin)), length(st))+1;

end
