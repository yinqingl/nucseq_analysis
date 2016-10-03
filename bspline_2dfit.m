%{
Point cloud fitting using bspline

Yinqing Li
yinqing@csail.mit.edu
2015

Minimize distance of point to s bspline and
minimize total length of the bspline or smoothness of the bspline

Usage:
[sp_fit,t_hat] = pc_bspline_fitting(x, sp_init, option)

Inputs:
x: point cloud, 2xn double array
sp_init: initial bspline
option: struct with fields
is_show_iteration = true;
is_show_f1_f2 = false;
norm_tol_th_iter = 1e-2;
n_delta_diff_f = 3;
max_iter = 300;
lambda = 1;
epsilon = 1e-10;
mode = d2' or 'tlen'
'd2' for minimizing secondary derivative
'tlen' for minimizing total length
    
Outputs:
sp_fit: bspline fitted to point cloud
t_hat: position of each point on the fitted bspline

Copyright (c) 2015, Yinqing Li
%}

function [sp_fit,t_hat] = bspline_2dfit(x, sp_init, option)

global epsilon
knots = sp_init.knots;
cp = sp_init.coefs;
cp = reshape(cp,numel(cp),1);
if size(x,1) > size(x,2),
    x = x';
end

if sp_init.order < 3,
    disp('sp_init.order should be larger than 3.');
    return
end

%normalize knots
s = max(max(x,[],2)-min(x,[],2)) / max(diff(knots));
knots = knots * s;

is_show_iteration = true;
is_show_f1_f2 = false;
norm_tol_th_iter = 1e-2;
n_delta_diff_f = 3;
max_iter = 300;
lambda = 1;
epsilon = 1e-10;
mode = 'tlen';
fixedse = true;

if nargin > 2,
    option_str = structvars(option);
    for i = [1:size(option_str,1)],
        eval(option_str(i,:));
    end
end

if strcmp(mode,'d2'),
    f = @f_d2;
    g = @g_d2;
elseif strcmp(mode,'tlen'),
    f = @f_tlen;
    g = @g_tlen;
end

if fixedse,
    Aeq = zeros(length(cp),length(cp));
    Aeq(1,1) = 1;
    Aeq(2,2) = 1;
    Aeq(end-1,end-1) = 1;
    Aeq(end,end) = 1;
    beq = zeros(size(cp));
    beq([1,2,end-1,end]) = cp([1,2,end-1,end]);
end
%iterate
t = knots(ceil(length(knots)/2)) * ones(size(x,2),1);
if_iter = true;
i_iter = 0;
diff_f_v = zeros(1,max_iter);
format shortg
tic

while(if_iter),
    cp_2d = reshape(cp,2,numel(cp)/2);
    sp = spmak(knots,cp_2d);
    
    t_hat = proj_x_to_sp(x,sp,t);
    f_cp = @(cp) f(x,knots,cp,t_hat,lambda);
    g_cp = @(cp) g(x,knots,cp,t_hat,lambda);
    f_g = @(cp) deal(f_cp(cp),g_cp(cp));

    %use fminunc
    %trust-region
    %quasi-newton
    if fixedse,
        cp_hat = fmincon(f_g,cp,[],[],Aeq,beq,[],[],[],optimoptions('fmincon','GradObj','on','Display','notify-detailed','UseParallel','always','Algorithm','interior-point'));
    else
        cp_hat = fminunc(f_g,cp,optimoptions('fminunc','GradObj','on','Display','notify-detailed','Algorithm','quasi-newton'));
    end
    
    %iteration condition
    i_iter = i_iter + 1;
    diff_f = norm(f_cp(cp_hat) - f_cp(cp));
    
    if i_iter > n_delta_diff_f,
        delta_diff_f = mean(abs(diff(diff_f_v(i_iter-n_delta_diff_f:i_iter-1))));
    else
        delta_diff_f = 0;
    end
    
    if diff_f < norm_tol_th_iter,  %convergence
        if_iter = false;
    elseif delta_diff_f < norm_tol_th_iter && i_iter > n_delta_diff_f,
        if_iter = false;
    elseif i_iter>=max_iter,
        if_iter = false;
    else   %continue iteration
        if is_show_iteration,
            disp(sprintf('iteration: %d',i_iter))
            disp(sprintf('diff_f = %f',diff_f))
            disp(sprintf('delta_diff_f = %f',delta_diff_f))
        end
        if is_show_f1_f2,
            [f_val,f1,f2] = f(x,knots,cp_hat,t_hat,lambda);
            disp(sprintf('f1 = %f, lambda*f2 = %f',f1,lambda*f2));
        end
        diff_f_v(i_iter) = diff_f;
    end
    cp = cp_hat;
    t = t_hat;
    toc
end

toc

knots = knots/s;
t_hat = t_hat/s;
sp_fit = spmak(knots,reshape(cp_hat,2,numel(cp_hat)/2));

end

function t_hat = proj_x_to_sp(x,sp,t_init)
t_hat = t_init;
for i = 1:size(x,2),
    y = @(t) pdist([x(:,i),fnval(sp,t)]','euclidean');
    t_hat(i) = fminsearch(y,t_init(i));
end
 
end

% guassian noise, smoothness based on secondary derivative
function [f_val,varargout] = f_d2(x,knots,cp,t_hat,lambda)
cp_2d = reshape(cp,2,numel(cp)/2);
sp = spmak(knots,cp_2d);
f1 = norm(x - fnval(sp,t_hat))^2;
f2 = norm(fnval(fnder( sp, 2),knots))^2;
f_val = f1 + lambda*f2;
if nargout > 1,
    varargout{1} = f1;
    varargout{2} = f2;
end
end

function g_val = g_d2(x,knots,cp,t_hat,lambda)
cp_2d = reshape(cp,2,numel(cp)/2);
sp = spmak(knots,cp_2d);
spd2 = fnder( sp, 2);
g_val_2d = zeros(size(cp_2d));
for i = 1:size(cp_2d,2),
    e_i = zeros(size(cp_2d));
    e_i(:,i) = 1;
    sp_e = spmak(knots,e_i);
    g1 = 2*sum((x - fnval(sp,t_hat)).*(-fnval(sp_e,t_hat)),2);
    spd2_e = fnder( sp_e, 2);
    g2 = 2*sum(fnval(spd2,knots).*fnval(spd2_e,knots),2);
    g_val_2d(:,i) = g1 + lambda*g2;
end
g_val = reshape(g_val_2d,numel(g_val_2d),1);
end

% guassian noise, smoothness based on total length
function [f_val,varargout] = f_tlen(x,knots,cp,t_hat,lambda)
cp_2d = reshape(cp,2,numel(cp)/2);
sp = spmak(knots,cp_2d);
f1 = norm(x - fnval(sp,t_hat))^2;

spd1 = fnder(sp,1);
d_sp = @(t) sqrt(sum(fnval(spd1,t).^2,1));
f2 = integral(d_sp,knots(1),knots(end));

f_val = f1 + lambda*f2;
if nargout > 1,
    varargout{1} = f1;
    varargout{2} = f2;
end
end

function g_val = g_tlen(x,knots,cp,t_hat,lambda)
global epsilon
cp_2d = reshape(cp,2,numel(cp)/2);
sp = spmak(knots,cp_2d);
g_val_2d = zeros(size(cp_2d));
for i = 1:size(cp_2d,2),
    e_i = zeros(size(cp_2d));
    e_i(:,i) = 1;
    sp_e = spmak(knots,e_i);
    g1 = 2*sum((x - fnval(sp,t_hat)).*(-fnval(sp_e,t_hat)),2);
    
    spd1 = fnder(sp,1);
    spd1_e = fnder( sp_e, 1);
    
    g2_int_xy = @(t) 1/2*repmat( (sum(fnval(spd1,t).^2,1) + epsilon).^(-1/2) ,2,1).*(2*fnval(spd1,t).*fnval(spd1_e,t));
    g2_int_x = @(t) cell2mat(cellfun(@(x) x(1,:), {g2_int_xy(t)}, 'UniformOutput', false));
    g2_int_y = @(t) cell2mat(cellfun(@(x) x(2,:), {g2_int_xy(t)}, 'UniformOutput', false));
%     g2 = integral(g2_int,knots(1),knots(end),'ArrayValued',true);
    g2 = [integral(g2_int_x,knots(1),knots(end));integral(g2_int_y,knots(1),knots(end))];
    
    g_val_2d(:,i) = g1 + lambda*g2;
end
g_val = reshape(g_val_2d,numel(g_val_2d),1);
end

% function x = g2_int_x(t,spd1,spd1_e)
% global epsilon
% g2 = 1/2*repmat( (sum(fnval(spd1,t).^2,1) + epsilon).^(-1/2) ,2,1).*(2*fnval(spd1,t).*fnval(spd1_e,t));
% x = g2(1,:);
% end
% 
% function y = g2_int_y(t,spd1,spd1_e)
% global epsilon
% g2 = 1/2*repmat( (sum(fnval(spd1,t).^2,1) + epsilon).^(-1/2) ,2,1).*(2*fnval(spd1,t).*fnval(spd1_e,t));
% y = g2(2,:);
% end

function example_code()
%use shortest distance to initialize bspline
label = double(dat_matrix_log_norm2(upper1lower('Slc12a5'),:));
figure; scatter(ydata(:,1),ydata(:,2),21,label(i_cells),'filled');
%use shortest path
d_th = 3;
ug_W = squareform(pdist(ydata,'euclidean'));
% k = 6+1;
for i = 1:size(ug_W,1),
%     [a,b] = sort(ug_W(i,:));
%     ug_W(i,b(k+1:end)) = 0;
    ug_W(i,ug_W(i,:)>d_th) = 0;
end
% ug_W_th = quantile(ug_W(:),1 - 6/(size(ydata,1)-1));
% ug_W(ug_W<ug_W_th) = 0;
ug_W = sparse(ug_W);

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
[sp_fit,t_hat] = pc_bspline_fitting(ydata,sp,option);
figure;
scatter(ydata(:,1),ydata(:,2),21,t_hat,'filled');
hold on
fnplt(sp_fit);
scatter(ydata_path(:,1),ydata_path(:,2),21,'filled');

%find distance between two ts along bspline
spd1 = fnder(sp_fit,1);
d_sp = @(t) sqrt(sum(fnval(spd1,t).^2,1));
d_t_sp = @(t1,t2) integral(d_sp,t1,t2);
t = linspace(0,1,100);
t = sort(t_hat);
d_t = zeros(size(t));
for i = 1:length(t),
    d_t(i) = d_t_sp(0,t(i));
end
figure;
plot(t,d_t,'.')

%show closest points on bspline
ydata_proj = fnval(sp_fit,t_hat);
figure;
scatter(ydata(:,1),ydata(:,2),21,t_hat,'filled');
hold on
fnplt(sp_fit);
scatter(ydata_path(:,1),ydata_path(:,2),21,'filled');
for i = 1:size(ydata,1),
    x = [ydata(i,1);ydata_proj(1,i)];
    y = [ydata(i,2);ydata_proj(2,i)];
    line(x,y,'color',[1 0 0],'LineWidth',2);
    hold on
end
daspect([1 1 1]);
end
