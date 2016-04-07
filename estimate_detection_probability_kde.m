%{
estimate the probability of a gene not expressed given that the gene is not observed

Yinqing Li
yinqingl@mit.edu
2015

inputs: 
dat: expression data matrix in tpm, 
     each row: genes, col: samples
option:
is_show_iteration: show details about minimization in each iteration
norm_tol_th_iter: iteration exiting condition, threshold for absolute changes in parameters
n_delta_diff_p_e_ij: the number of iterations whose parameter changes are averaged
max_iter: maximum iteration
n_kde: the number of bins of expression levels for estimating expression distribution 
epsilon: a small number for numerical stability in division
beta_lb: lower bound for parameters used in geometric curve fitting
alg: algorithm, 'geometric','geometric_log','logit','bspline'

output: 
p_det: associated probability estimation
use parametric model to fit detection rate for cell i
use kde to fit distribution of gene j
p_det is returned as sparse matrix
p_det = ones(size(dat))
p_det(dat==0) = p_det;

Copyright (c) 2015, Yinqing Li
%}

function [p_det,varargout] = estimate_detection_probability_kde(dat,option)

if max(max(dat)) < 100,
    disp('dat should be on linear scale.')
    p_det = [];
    varargout = {};
    return
end

% tic
%dat format, rows: genes, columns: samples
%dat should not contain zero row

is_show_iteration = true;
norm_tol_th_iter = 1e-2;
n_delta_diff_p_e_ij = 3;
max_iter = 300;
n_kde = 2^8;
epsilon = 1e-4;
beta_lb = [-epsilon,epsilon];    %geometric
alg = 'geometric'; %{'geometric','geometric_log','logit','bspline'}

if nargin >1,
    option_str = structvars(option);
    for i = [1:size(option_str,1)],
        eval(option_str(i,:));
    end
end

%disp options
disp(sprintf('show iteration: %d',is_show_iteration))
disp(sprintf('algorithm: %s',alg))
disp(sprintf('number of kde bins: %d',n_kde))

%remove zero rows
j_z = sum(dat~=0,2)==0;
dat_in = dat;
dat = dat(~j_z,:);

% load tmp_dat
[n_j, n_i] = size(dat);

i_nz = (dat~=0);
i_z = (dat == 0);

%initial condition
epsilon = 1e-3;
p_e_ij_z_d_ij_z = ones(size(dat));
p_e_ij_z_d_ij_z(i_z) = 0.5;  %0.5 probability e  = 0 when see d = 0;

if isfield(option, 'prev_p_det'),
    p_e_ij_z_d_ij_z(:,i_X_ref) = prev_p_det(~j_z,:);
    disp('use prev_p_det');
end

p_d_ij_nz_e_ij_nz = ones(size(dat)) - epsilon;   %prob d_ij ~= 0 when e_ij ~= 0
p_d_ij_z_e_ij_nz = 1-p_d_ij_nz_e_ij_nz;   %prob d_ij = 0 when e_ij ~= 0

%estimate p_e_j for gene j for its nonzero kde distribution
%distribution is estimated on log scale
dat_log = log(dat+1);
p_e_j = struct();
% p_e_j.e = zeros(n_j,n_kde+1);  %support
p_e_j.e = zeros(1,n_kde+1);  %support
p_e_j.p = zeros(n_j,n_kde+1);  %distribution

%initialize with naive bayes
p_e_j_z = sum(p_e_ij_z_d_ij_z.*i_z,2)./sum(p_e_ij_z_d_ij_z,2);
p_e_j_z(p_e_j_z < epsilon) = epsilon;
p_e_j_z(p_e_j_z > 1 - epsilon) = 1 - epsilon;

%use the same support for all genes to estimate distribution
minimum=min(dat_log(i_nz)); maximum=max(dat_log(i_nz));
Range=maximum-minimum;
MIN=minimum-Range/10; MAX=maximum+Range/10;

tic
display('start kernel density estimation');
for j = [1:n_j],
    
    dat_log_j_nz = dat_log(j,i_nz(j,:));
    %{
    [~,density,xmesh,~]=kde(dat_log_j_nz,n_kde);
    density = density/sum(density);
    %}
    
    [~,density,xmesh,~]=kde(dat_log_j_nz,n_kde,MIN,MAX);
    density(density<0) = 0;
    density = density/sum(density);
    
    i_z_xmesh = (xmesh<=0);
    xmesh(i_z_xmesh) = 0;
    density(i_z_xmesh) = 0;    
    
%     p_e_j.e(j,:) = [0,xmesh];
    p_e_j.p(j,:) = [p_e_j_z(j),(1-p_e_j_z(j))*density'];
end
toc
% p_e_j.e = exp(p_e_j.e)-1;
p_e_j.e = exp([0,xmesh])-1;
p_e_j.e = reshape(p_e_j.e,numel(p_e_j.e),1);

%lambda_i cannot be initialized with zero
epsilon = 1e-4;

%logit
if strcmp(alg,'logit'),
    beta_i = repmat([-5;1],1,n_i);  %logit
    lambda = @lambda_logit;
    d_lambda = @d_lambda_logit;
    data = dat_log; %fit on log scale
    p_e_j.e = log(p_e_j.e+1); %fit on log scale
    beta_A = [];
    beta_b = [];    
    beta_lb = [-inf,epsilon]; %logit
    beta_ub = [-epsilon,inf]; %logit
end

%geometric
if strcmp(alg,'geometric'),
    e_s = 10^(-round(log10(max(max(dat_in)))/2));
    disp(sprintf('for numerical stability, scale e by %.1e.',e_s));
    
    beta_i = repmat([0;2.5],1,n_i);   %geometric
    lambda = @(e,beta) lambda_geometric(e*e_s,beta); %numerical stability
    d_lambda = @(e,beta) d_lambda_geometric(e*e_s,beta);
    data = dat; %fit on linear scale
    %p_e_j.e is on linear scale
    beta_A = [];
    beta_b = [];
    % beta_lb = [-epsilon,epsilon];    %geometric
    beta_ub = [0,inf];    %geometric
end

%geometric_log
if strcmp(alg,'geometric_log'),
    beta_i = repmat([0;0.25],1,n_i);   %geometric
    lambda = @lambda_geometric;
    d_lambda = @d_lambda_geometric;
    data = dat_log; %fit on log scale
    p_e_j.e = log(p_e_j.e+1); %fit on log scale
    beta_A = [];
    beta_b = [];
    % beta_lb = [-epsilon,epsilon];    %geometric
    beta_ub = [0,inf];    %geometric
end

%bspline
if strcmp(alg,'bspline'),
    %initialize bspline for each cell
    n_cp = 12;
    cp = linspace(0,1,n_cp)';
    % cp = [0,0.01,0.02,0.5,0.9,0.95,0.96]';
    k = 3+1;
    dat_log_max = max(dat_log,[],1) + 1;
    beta_i = repmat(cp,1,n_i);
    m_knots = zeros(n_cp+k,n_i);
    for i = 1:n_i,
        m_knots(:,i) = augknt(linspace(0,dat_log_max(i),n_cp-2)',k);
    end
    data = dat_log; %fit on log scale
    p_e_j.e = log(p_e_j.e+1); %fit on log scale
    beta_A = zeros(n_cp-1,n_cp);
    beta_A(1:n_cp:end) = 1;
    beta_A(n_cp:n_cp:end) = -1;
    beta_b = zeros(n_cp-1,1);
    beta_lb = zeros(n_cp,1);
    beta_ub = ones(n_cp,1);
    beta_ub(1) = 1e-5;  %the detection prob for 0 is 0
end

if isfield(option, 'prev_beta_i'),
    beta_i(:,i_X_ref) = prev_beta_i;
    disp('use prev_beta_i');
end
if isfield(option, 'prev_knots_i'),
    m_knots(:,i_X_ref) = prev_knots_i;
    disp('use prev_knots_i');
end

% lambda = lambda_geometric;
% lambda = @lambda_logit;
% lambda = @lambda_geometric;
% @lambda_logit(e,beta)
% @lambda_geometric(e,beta)

%iterate
if_iter = true;
i_iter = 0;
diff_p_e_ij_v = [];
format shortg
tic
while(if_iter),
    
    %estimate lambda_i
    for i = [1:n_i],
%         tic       
        e_ij = data(:,i);
        if strcmp(alg,'bspline'),
            lambda = @(e,beta) lambda_bspline(e,beta',m_knots(:,i));
            d_lambda = @(e,beta) d_lambda_bspline(e,beta',m_knots(:,i));
        end           
        f = @(x) -f_ll_numel(e_ij, p_e_j, x, lambda);       
        g = @(x) -g_ll(e_ij, p_e_j, x, lambda, d_lambda);
        
        f_g = @(x) deal(f(x),g(x));

    %     %use fzero
    %     lambda_ik(k,i) = fzero(g,lambda_ikk(k,i));
    % 
    %     %use fminsearch
    %     lambda_ik(k,i) = fminsearch(f,lambda_ik_est(k,i),optimset('MaxFunEvals',1e3,'Display','off'));

        %use fmincon
        %trust-region-reflective
        %interior-point
        beta_i(:,i) = fmincon(f_g,beta_i(:,i),beta_A,beta_b,[],[],beta_lb,beta_ub,[],optimoptions('fmincon','GradObj','on','Display','notify-detailed','UseParallel','always','Algorithm','interior-point'));
%         beta_i(:,i) = fmincon(f_g,beta_i(:,i),[],[],[],[],beta_lb,beta_ub,[],optimoptions('fmincon','GradObj','on','Display','notify-detailed','UseParallel','always','Algorithm','trust-region-reflective'));
%         if rem(i,10) == 0,
%             i
%             c = strrep(num2str(clock),' ','');
%             disp(c)
%             toc
%         end
    end
    
    toc
    
    %detection probability
    s_e_ij = p_e_j.e(2:end); %support for distribution
    
    if strcmp(alg,'bspline'),
        for i = [1:n_i],
            lambda = @(e,beta) lambda_bspline(e,beta',m_knots(:,i));
            p_d_ij_z_e_ij_nz(:,i) = p_e_j.p(:,2:end)*(1-lambda(s_e_ij,beta_i(:,i)));
        end
    else
        p_d_ij_z_e_ij_nz = p_e_j.p(:,2:end)*(1-lambda(s_e_ij,beta_i));
    end
    
%     for i =[1:n_i]
%         p_d_ij_z_e_ij_nz(i_z(:,i),i) = p_e_j.p(i_z(:,i),2:end)*(1-lambda(s_e_ij,beta_i(:,i))); %prob d_ij = 0 when e_ij ~= 0
%     end

    %estimate p_e_ij_d_ij
    p_e_ij_z = repmat(p_e_j.p(:,1),1,n_i);
    p_e_ij_z_d_ij_z(i_z) = p_e_ij_z(i_z)./(p_e_ij_z(i_z) + p_d_ij_z_e_ij_nz(i_z));
    p_e_ij_z_d_ij_z(i_nz) = 1;
    
    
    %update p_e_j
    epsilon = 1e-3;
    p_e_j_z = sum(p_e_ij_z_d_ij_z.*i_z,2)./sum(p_e_ij_z_d_ij_z,2);
    p_e_j_z(p_e_j_z < epsilon) = epsilon;
    p_e_j_z(p_e_j_z > 1 - epsilon) = 1 - epsilon;
    
    p_e_j_p_previter = p_e_j.p;
    ita = (1-p_e_j_z)./(1-p_e_j_p_previter(:,1));
    for j = [1:n_j]
        p_e_j.p(j,:) = [p_e_j_z(j),ita(j)*p_e_j_p_previter(j,2:end)];
    end
    
    %iteration condition
    i_iter = i_iter + 1;
    diff_p_e_ij = norm(p_e_j_p_previter(:,1) - p_e_j.p(:,1));
    
    if length(diff_p_e_ij_v) > n_delta_diff_p_e_ij,
        delta_diff_p_e_ij = mean(abs(diff(diff_p_e_ij_v(end-n_delta_diff_p_e_ij+1:end))));
    else
        delta_diff_p_e_ij = 0;
    end
    
    if diff_p_e_ij < norm_tol_th_iter,  %convergence
        if_iter = false;
    elseif delta_diff_p_e_ij < norm_tol_th_iter && length(diff_p_e_ij_v) > n_delta_diff_p_e_ij,
        if_iter = false;
    elseif i_iter>max_iter,
        if_iter = false;
    else,   %continue iteration
        if is_show_iteration,
            disp(sprintf('iteration: %d',i_iter))
            disp(sprintf('diff_p_e_ij = %f',diff_p_e_ij))
            disp(sprintf('delta_diff_p_e_ij = %f',delta_diff_p_e_ij))
        end
        diff_p_e_ij_v(end+1) = diff_p_e_ij;
    end
    toc
end

%p_det is defined as the probability of e_ij == d_ij
% p_det = zeros(size(dat));
% p_det(i_nz) = 1;
% p_det(i_z) = p_e_ij_z_d_ij_z(i_z);
p_det = ones(size(dat_in));
p_det(~j_z,:) = p_e_ij_z_d_ij_z;
i_z = dat_in == 0;
p_det = p_det(i_z);

%show iteration
%{
format shortg
c = clock;
c = strrep(num2str(fix(c)),' ','');

figure();
plot(diff_p_e_ij_v,'-.');
plotname = sprintf('lambda_iteration_iter=%d_%s',max_iter,c);
title(plotname,'Interpreter','none')
xlabel('iteration')
ylabel('diff_p_e_ij')
print(gcf,'-r250','-dpng',['analysis_figures/' plotname '.png'])

%show lambda
figure();
cmap = lines(n_i+5);
e = exp(linspace(min(dat_log(:)),max(dat_log(:)),100))-1;
for i = [1:n_i],
    plot(log(e+1),lambda(e',beta_i(:,i)),'-.','color',cmap(i,:))
    hold on
end
plotname = sprintf('lambda_iter=%d_%s',max_iter,c);
title(plotname,'Interpreter','none')
xlabel('log(e_ij+1)','Interpreter','none')
ylabel('lambda_i','Interpreter','none')
print(gcf,'-r250','-dpng',['analysis_figures/' plotname '.png'])
%}
prob_est = struct();
if strcmp(alg,'bspline'),
	lambda = @lambda_bspline;
    prob_est.m_knots = m_knots;
end
prob_est.lambda = lambda;
prob_est.beta_i = beta_i;
prob_est.p_e_j = p_e_j;
varargout{1} = prob_est;

% figure()
% imagesc(log(dat(find(sum(dat,2)),:))+1);
% figure()
% imagesc(p_det(find(sum(dat,2)),:));
end

%log likelihood function
%probabilty for all genes defined on the same support (min(tpm), max(tpm))
%plug in lambda as a function
function ll = f_ll_numel(e_ij, p_e_j, beta_i, lambda)
i_nz = e_ij>0;
i_z = e_ij == 0;

% s_e_ij = p_e_j.e(i_z,2:end); %support for distribution of a gene
s_e_ij = p_e_j.e(2:end); %support for distribution
p_z = p_e_j.p(i_z,2:end);    %the first support is zero
p_nz = 1-p_e_j.p(i_nz,1);

% ll = sum(log(lambda(e_ij(i_nz),beta_i))) + sum(log(1-sum(p.*lambda(s_e_ij,beta_i),2)));
ll = sum( log( lambda(e_ij(i_nz),beta_i) + realmin ) ) + ...
    sum( p_z *  log(1-lambda(s_e_ij,beta_i) + realmin) );
%     sum( log(p_nz) ) + ...    %constant with respect to beta
%     sum( p_z(:) .* log(p_z(:)+realmin) );
end

% function d_ll = g_ll_numel(e_ij, p_e_j, beta_i, lambda)
% e = 1e-5;
% e_v = e*eye(length(beta_i));
% d_ll = size(beta_i);
% for i=[1:size(e_v,2)],
%     d_ll(i) = (f_ll_numel(e_ij, p_e_j, beta_i+e_v(:,i), lambda)-f_ll_numel(e_ij, p_e_j, beta_i-e_v(:,i), lambda))/2/e;
% end
% end

function d_ll = g_ll(e_ij, p_e_j, beta_i, lambda, d_lambda)
i_nz = e_ij>0;
i_z = e_ij == 0;

% s_e_ij = p_e_j.e(i_z,2:end); %support for distribution of a gene
s_e_ij = p_e_j.e(2:end); %support for distribution
p_z = p_e_j.p(i_z,2:end);    %the first support is zero

d_lambda_e_ij_nz = d_lambda(e_ij(i_nz),beta_i);
lambda_e_ij_nz = lambda(e_ij(i_nz),beta_i) + 1e-6;   %add 1e-5
d_lambda_e_ij_z = d_lambda(s_e_ij,beta_i);
lambda_e_ij_z = lambda(s_e_ij,beta_i);

d_ll = size(beta_i);
% d_ll(1) = sum(d_lambda_e_ij_nz./lambda_e_ij_nz) + ...
%     sum( p_z * (-d_lambda_e_ij_z./(1-lambda_e_ij_z + 1e-6)) );
% d_ll(2) = sum(d_lambda_e_ij_nz.*e_ij(i_nz)./lambda_e_ij_nz) + ...
%     sum( p_z * (-s_e_ij.*d_lambda_e_ij_z./(1-lambda_e_ij_z + 1e-6)) );
for i = 1:length(beta_i),
    d_ll(i) = sum(d_lambda_e_ij_nz(:,i)./lambda_e_ij_nz) + ...
        sum( p_z * (-d_lambda_e_ij_z(:,i)./(1-lambda_e_ij_z + 1e-6)) );
end
end

function lambda = lambda_logit(e,beta)
% [n_row, n_col] = size(e);
% e = reshape(e,1,numel(e));
% e = log(e+1); %fit on log scale
lambda = 2 * (1./(1+exp(-[ones(size(e)),e]*beta)) - 0.5);
% lambda = reshape(lambda,n_row,n_col);
end

function d_lambda = d_lambda_logit(e,beta)
% [n_row, n_col] = size(e);
% e = reshape(e,1,numel(e));
% e = log(e+1); %fit on log scale
t = [ones(size(e)),e]*beta;
exp_t = exp(-t);
dlambda_dt = 2*exp_t./(1+exp_t).^2;

dt_dbeta = [ones(size(e)),e];
d_lambda = zeros(length(e),length(beta));
for i = 1:length(beta),
    d_lambda(:,i) = dlambda_dt.*dt_dbeta(:,i);
end

% d_lambda = reshape(d_lambda,n_row,n_col);
end

%beta > 0
function lambda = lambda_geometric(e,beta)
% [n_row, n_col] = size(e);
% e = reshape(e,1,numel(e));
% e = log(e+1); %fit on log scale

t = [ones(size(e)),e]*beta;
x = 1-exp(-t);
lambda = max(x,0);  %rectifier

% lambda(i) = log(exp(x(i)*2e1)+1)/2e1;
% lambda = log(exp(x*2e1)+1)/2e1;
% lambda = 1-exp(-t);
% lambda = reshape(lambda,n_row,n_col);
end

function d_lambda = d_lambda_geometric(e,beta)
% [n_row, n_col] = size(e);
% e = reshape(e,1,numel(e));
% e = log(e+1); %fit on log scale
t = [ones(size(e)),e]*beta;

% s_t = t<0;
% d_lambda = zeros(size(t));
% d_lambda(~s_t) = exp(-t(~s_t));
% y=1/4*(2+e.^x-(e.^x-2*x).*sech(x).^2-(-2+e.^x).*tanh(x));

x = 1-exp(-t);
dlambda_dt = 1./(1+exp(-1e2*x)).*exp(-t);
dt_dbeta = [ones(size(e)),e];
d_lambda = zeros(length(e),length(beta));
for i = 1:length(beta),
    d_lambda(:,i) = dlambda_dt.*dt_dbeta(:,i);
end
% d_lambda = exp(-t);

% d_lambda = reshape(d_lambda,n_row,n_col);
end

function lambda = lambda_bspline(e,beta,knots)
cp = beta;
sp = spmak(knots,cp);
lambda = fnval(sp,e);
end

function d_lambda = d_lambda_bspline(e,beta,knots)
cp = beta;
d_lambda = zeros(length(e),length(cp));
for i = 1:length(cp),
    e_i = zeros(size(cp));
    e_i(:,i) = 1;
    sp_e = spmak(knots,e_i);
    d_lambda(:,i) = fnval(sp_e,e);
end

end

function test_function()
import bioma.data.*

%import functions
fh = custom_analysis_functions;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

load dat_20150414_norm_rDESeq2_DEseq2

dat_matrix_norm = dat_matrix_norm2;
dat_matrix_log_norm = log(dat_matrix_norm+1);

[dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(dat_matrix_norm);

fh = functions_test;
fh_str = structvars(fh);
for i = [1:size(fh_str,1)],
    eval(fh_str(i,:));
end

option = struct();
option.is_show_iteration = true;
option.max_iter = 300;

%logit
x = exp(linspace(0,MAX,100))-1;
figure;plot(log(x+1),lambda(log(x+1)',beta_i(:,i)));

%geometric
x = exp(linspace(0,MAX,100))-1;
figure;plot(log(x+1),lambda(x',beta_i(:,i)));

%visualize bspline
figure;
for i = 1:n_i,
    x = linspace(0,max(m_knots(:,i)),100);
    plot(x,lambda_bspline(x,beta_i(:,i)',m_knots(:,i)),'LineWidth',2); hold on;
end


end


%{
%log likelihood function
%probabilty for all genes defined on the same support (min(tpm), max(tpm))
%e_ij should be rounded to the closet point on the finite support defined on log scale
%I_e_ij is the index of transformed e_ij on the finite support
function ll = f_ll(I_e_ij, p_e_j, lambda_e_i_beta_i)
lambda_e_i = lambda_e_i_beta_i(1:end-2);
beta_i = lambda_e_i_beta_i(end-1:end);  %first order polynomial
i_nz = e_ij>0;
i_z = e_ij == 0;

lambda_e_ij = lambda_e_i(I_e_ij);

p = p_e_j.p(:,2:end);    %the first support is zero

ll = sum(log(lambda_e_ij(i_nz))) + sum(log(1-p(i_z)*lambda_e_ij));
end

%gradient of log likelihood function
function d_ll = g_ll(I_e_ij, p_e_j, lambda_e_i_beta_i)
lambda_e_i = lambda_e_i_beta_i(1:end-2);
beta_i = lambda_e_i_beta_i(end-1:end);  %first order polynomial
i_nz = e_ij>0;
i_z = e_ij == 0;

lambda_e_ij = lambda_e_i(I_e_ij);

p = p_e_j.p(:,2:end);    %the first support is zero

d_ll = sum(1./lambda_e_ij(i_nz)) - sum(p,2)./(1-p(i_z)*lambda_e_ij);
end

function d_ll_numel = g(I_e_ij, p_e_j, lambda_e_i_beta_i)
e = 1e-5;
e_v = e*eye(length(lambda_e_i_beta_i));
for i=[1:size(e_v,2)],
    d_ll_numel = (f(I_e_ij, p_e_j, lambda_e_i_beta_i+e_v(:,i)) - f(I_e_ij, p_e_j, lambda_e_i_beta_i-e_v(:,i)))/2/e;
end
end

%nonlinear constraints
%logistic regression
function [c,ceq,d_c,d_ceq] = nlcon_logit(I_e_ij, p_e_j, lambda_e_i_beta_i)
lambda_e_i = lambda_e_i_beta_i(1:end-2);
beta_i = lambda_e_i_beta_i(end-1:end);  %first order polynomial
i_nz = e_ij>0;
i_z = e_ij == 0;

lambda_e_ij = lambda_e_i(I_e_ij);

e_ij = p_e_j.e(I_e_ij); %support for distribution of a gene

c = [];
ceq = sum((log(lambda_e_ij./(1-lambda_e_ij)) - beta_i'*[ones(size(e_ij));e_ij]).^2);
d_c = [];
d_ceq = g_nlcon_logit(I_e_ij, p_e_j, lambda_e_i_beta_i);
end

function d_ceq = g_nlcon_logit(I_e_ij, p_e_j, lambda_e_i_beta_i)
lambda_e_i = lambda_e_i_beta_i(1:end-2);
beta_i = lambda_e_i_beta_i(end-1:end);  %first order polynomial
i_nz = e_ij>0;
i_z = e_ij == 0;

lambda_e_ij = lambda_e_i(I_e_ij);

e_ij = p_e_j.e(I_e_ij); %support for distribution of a gene

c_ij = (log(lambda_e_ij./(1-lambda_e_ij)) - beta_i'*[ones(size(e_ij));e_ij]);

d_lambda_e_ij = 2*sum(c_ij.*(1./lambda_e_ij./(1-lambda_e_ij)-beta_i(0)));
d_beta_0 = -2*sum(c_ij);
d_beta_1 = -2*sum(c_ij.*e_ij);
d_ceq = [d_lambda_e_ij;d_beta_0;d_beta_1];
end
%}