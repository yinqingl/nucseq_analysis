%{
help functions for visualizing the selected bisne genes

Yinqing Li
yinqingl@mit.edu
MIT, 2015
%}

function scatter_series(ydata, dat_matrix)
figure;
global ii
ii = 1;
if ~strcmp(class(dat_matrix),'bioma.data.DataMatrix')
    import bioma.data.*
    dat_matrix = double(dat_matrix);
    dat_matrix = DataMatrix(dat_matrix, 1:size(dat_matrix,1), 1:size(dat_matrix,2));
end
    
[dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(dat_matrix);
scatter(ydata(:,1),ydata(:,2),21,dat(ii,:),'filled');
title(sprintf('%d %s',ii,genes{ii}));
dat_in = {ydata, dat_matrix};
set(gcf, 'WindowButtonDownFcn', {@mouse_click_update_ii, dat_in});
set(gcf,'KeyPressFcn', {@key_press_update_ii, dat_in});
end

function mouse_click_update_ii(object, eventdata, dat_in)
%left click
% persistent ii;
% if isempty(ii),
%     ii = 1;
% end
global ii
ydata = dat_in{1};
dat_matrix = dat_in{2};
if strcmp(get(gcbf, 'SelectionType'),'normal'),
    cp = get(gca, 'CurrentPoint');
    x = round(cp(1,1));
%     y = round(cp(1,2));
    if x < median(ydata(:,1)),
        ii = ii - 1;
        if ii <= 0
            ii = size(dat_matrix,1);
        end
    else
        ii = ii + 1;
        if ii > size(dat_matrix,1);
            ii = 1;
        end
    end
    [dat, cellnames, genes, num_samples, num_genes] = datamatrix2matrix(dat_matrix);
    scatter(ydata(:,1),ydata(:,2),21,dat(ii,:),'filled');
    title(sprintf('%d %s',ii,genes{ii}));
    drawnow
end

end