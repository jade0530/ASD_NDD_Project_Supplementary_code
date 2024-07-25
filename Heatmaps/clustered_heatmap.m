
%%%%%%%%%%% read the matrix %%%%%%%%%%%
inputfile='../inputdata.txt';
fid=fopen(inputfile, 'r');
formatstring='%s %s %s %s %s %s %s %s %s %s %s %s %s %*[^\n]';
heatmap_matrix=textscan(fid, formatstring, 'delimiter', '\t');

%%%%%%%%%%% generate clustered heatmap %%%%%%%%%%%
cg = clustergram(heatmap_matrix,'RowLabels', xvalues_ASD, 'ColumnLabels', yvalues, 'Colormap',redbluecmap, 'Standardize','row');

%%%%%%%%%%% Get all handles from root %%%%%%%%%%%
set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');

%%%%%%%%%%% Find hearmap axis and change the font size %%%%%%%%%%%
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 6)