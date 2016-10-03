%clustering based on local density
%{
Clustering by fast search and find of density peaks
Alex Rodriguez, Alessandro Laio

modified by
Yinqing Li
yinqingl@csail.mit.edu
11/4/2014

cl: cluster assignment
stat: structure
stat.cluster : number of clusters
stat.center : center of clusters
stat.elements : number of elements in each cluster
stat.core : number of core elements in each cluster, defined by percent
%}

function [cl,stat] = cluster_dp_pdist(dat)
xx=pdist(dat, 'euclidean');
ND = size(dat,1);
N=length(xx);
dist = squareform(xx);

percent=2.0;
fprintf('Average percentage of neighbours (hard coded): %5.6f\n', percent);

position=round(N*percent/100);
sda=sort(xx);
dc=sda(position);

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);

rho = zeros(ND,1);
%
% Gaussian kernel
%
rho = sum(exp(-(dist/dc).^2),2)-1;
%{
for i=1:ND-1
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)^2);
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end
%}

%
% Average Gaussian kernel
%
rho_n = sum(dist<2*dc,2)-1;
rho = rho./rho_n;

%
% "Cut off" kernel
%
% rho = sum(dist<dc,2)-1;
%{
for i=1:ND-1
 for j=i+1:ND
   if (dist(i,j)<dc)
      rho(i)=rho(i)+1.;
      rho(j)=rho(j)+1.;
   end
 end
end
%}

maxd=max(dist(:));

[rho_sorted,ordrho]=sort(rho,'descend');
delta = zeros(length(ordrho),1);
nneight = zeros(length(ordrho),1);
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
ind = [1:ND];
gamma = rho.*delta;
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')

rect = getrect;
rhomin=rect(1);
deltamin=rect(2);
NCLUST=0;
cl = ones(ND,1)*(-1);
for i=1:ND
  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;
     icl(NCLUST)=i;
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end

%halo, core of each cluster
%halo ~= 0, indicates in core
for i=1:ND
  halo(i)=cl(i);
end
if (NCLUST>1)
  for i=1:NCLUST
    bord_rho(i)=0.;
  end
  for i=1:ND-1
    for j=i+1:ND
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end
nc_cellarray = {};
nh_cellarray = {};
for i=1:NCLUST
  nc=0;
  nh=0;
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
  nc_cellarray{i} = nc;
  nh_cellarray{i} = nh;
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end
stat = struct();
stat.cluster = i;
stat.center = icl;
stat.elements = nc_cellarray;
stat.core = nh_cellarray;

%{
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
disp('Performing 2D nonclassical multidimensional scaling')
Y1 = mdscale(dist, 2, 'criterion','metricstress');
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (cl(j)==i) %only plot the data in core of each cluster
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
%}

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end

%{
faa = fopen('CLUSTER_ASSIGNATION', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
%}