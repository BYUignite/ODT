
path(path,'../tools');

%----------------- parameters

nbins   = 40;                             % number of uniform bins in PDF
ntpdf   = 40;                             % number of time windows for PDF data

%------------------ Initialize variables

command = '../tools/getInputFileParameter.py ../../input/odtParam.inp tEnd';
[status tEnd] = system(command);
tEnd = str2num(tEnd);

Py      = zeros(nbins,ntpdf);       % var PDF bined in nbins of var, for number of time windows, and npart
Pl      = zeros(nbins,ntpdf);       % var PDF bined in nbins of var, for number of time windows, and npart
tpdfw   = linspace(0,tEnd,ntpdf+1);       % tpdfw(1)-->tpdfw(2) is window 1; tpdf(2)-->tpdfw(3) is window 3 etc
tpdf    = (tpdfw(2:end)+tpdfw(1:end-1))/2;% times of the pdf window centers

%------------------  Read in data

fileName = 'eddyData.dat';
data = readMatData(fileName);

t_eddy = data(:,1);
y_eddy = data(:,2);
l_eddy = log10(data(:,3));

llo=min(l_eddy);
lhi=max(l_eddy);

ylo=min(y_eddy);
yhi=max(y_eddy);


%------------------ Make the pdf

for itw = 1:ntpdf                     % loop pdf time windows
    ii = find(t_eddy >= tpdfw(itw) & t_eddy <= tpdfw(itw+1));  % time indicies of t_eddy in pdf time window
    xx1 = y_eddy(ii);         % all data points
    xx2 = l_eddy(ii);         % all data points
    [ybins Py(:,itw)] = compute_pdf(xx1,nbins,ylo,yhi);
    [lbins Pl(:,itw)] = compute_pdf(xx2,nbins,llo,lhi);
end

%------------------ Write the files

fileName = strcat('Peddy_y.dat');
fid = fopen(fileName, 'w');
fprintf(fid, '# times -->\n# ');
fprintf(fid, '%f  ', tpdf);
fprintf(fid, '\n# y_(m), Pdf...\n');
data = [ybins Py];
[ni nj] = size(data);
for i=1:ni
    fprintf(fid, '%-16.8e', data(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

fileName = strcat('Peddy_l.dat');
fid = fopen(fileName, 'w');
fprintf(fid, '# times -->\n# ');
fprintf(fid, '%f  ', tpdf);
fprintf(fid, '\n# l_(m), Pdf...\n');
data = [lbins Pl];
[ni nj] = size(data);
for i=1:ni
    fprintf(fid, '%-16.8e', data(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

%------------------ Plot data

hFig = figure('visible', 'off');
clf;

subplot(1,2,1);
[T,Y]=meshgrid(tpdf,ybins);
contourf(Y,T,Py,50);
shading flat;
title('Eddy Location PDF', 'FontSize', 16);
xlabel('Position (m)', 'FontSize', 16);
ylabel('Time (s)', 'FontSize', 16);
set(gca,'FontSize',16);

subplot(1,2,2);
[T,L]=meshgrid(tpdf,lbins);
contourf(L,T,Pl,50);
shading flat;
title('Eddy Size PDF', 'FontSize', 16);
xlabel('log_{10}(Eddy Size/(m))', 'FontSize', 16);
ylabel('Time (s)', 'FontSize', 16);
set(gca,'FontSize',16);
%h=colorbar('East','FontSize',16);
%set(h,'YColor',[1 1 1]);

hgexport(gcf, 'eddySizePosPdfContours.pdf', hgexport('factorystyle'), 'Format', 'pdf');

exit;

