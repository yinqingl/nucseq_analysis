%{
smooth using lowess method and fit data with glm

Yinqing Li
yinqing@mit.edu
2015

input:
x: vector
y: vector
span: between 0 and 1, the portion of x used in local regression
iter: number of iterations

output:
yest: estimated y

Copyright (c) 2015, Yinqing Li
%}

function yest = lowess_glm(x, y, span, iter)

if nargin < 4,
    iter = 3;
end

i_nan = isnan(x) | isnan(y) | (y<=0);
n_nan = length(x);

x = x(~i_nan);
y = y(~i_nan);

n = length(x);
r = ceil(span*n);
h = zeros(1,n);
tic
for i = [1:n],
   tmp = sort(abs(x - x(i))); %sort ascend by default
   h(i) = tmp(r)+1e-9;
end
toc
yest = zeros(size(y));
delta = ones(size(x));
for iteration = [1:iter],
    for i = [1:n],
        %for memory reason, calculate w at each iteration
        %if matrix is small, calculate and store w before iteration
        w = clip(abs(x-x(i))/h(i),0.0,1.0);
        w = (1-w.^3).^3;
        
        %adjust span if data is too much
        %remove ill-conditioned weights
        max_num_dat_glm = 2000;
        w_th = 1e-2;
        w_index = (w>w_th);
        if sum(w_index) > max_num_dat_glm,
            [~,i_s_w] = sort(w(w_index));
            w_index = w_index(i_s_w(end-max_num_dat_glm+1:end));
        end
        
        w = w(w_index);
        weights = delta(w_index).*w;
        
        %perform glm with weight
        warning off
        try,
            b = glmfit(x(w_index),y(w_index),'gamma','weights',weights);
        catch ME,
            1;
        end
        yest(i) = glmval(b,x(i),'reciprocal');
    end
    residuals = y - yest;
    s = median(abs(residuals));
    delta = clip(residuals/(6*s),-1,1);
    delta = (1-delta.^2).^2;
    toc
end

yest_nan = nan(1,n_nan);
yest_nan(~i_nan) = yest;
yest = yest_nan;
end

function x = clip(x,min,max)
x(x<min)=min;
x(x>max)=max;
end

