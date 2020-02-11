
% nDmp is passed in as an argument

nbins = 50;

path(path,'../tools');

Pdx = zeros(nbins,nDmp);

command = '../tools/getInputFileParameter.py ../../input/odtParam.inp domainLength';
[status Ld] = system(command);
Ld = str2num(Ld);

command = '../tools/getInputFileParameter.py ../../input/odtParam.inp dxmin';
[status dxmin] = system(command);
dxmin = str2num(dxmin) * Ld;

command = '../tools/getInputFileParameter.py ../../input/odtParam.inp dxmax';
[status dxmax] = system(command);
dxmax = str2num(dxmax) * Ld;

dxmin = log10(dxmin);
dxmax = log10(dxmax);

xmin = dxmin-0.5;
xmax = dxmax+0.5;

%----------------

for idmp = 1:nDmp

    idmp

    fname = strcat('gridSizes/gridSizes_', num2str(idmp), '.dat');
    data = readMatData(fname);

    xf = data(:,1);
    dx = xf(2:end)-xf(1:end-1);
    ifix = find(dx < 0.0);
    dx(ifix) = Ld - xf(ifix);

    dx = log10(dx);

    [xbins Pdx(:,idmp)] = compute_pdf(dx,nbins,xmin,xmax);

end

%---------------- grab the dump times

dumpTimes = [];

fname = '../../input/dumpTimes.inp';
fid = fopen(fname, 'r');
ln = strtrim( fgetl(fid) );
i = 1;
while(~feof(fid))
    ln = strtrim( fgetl(fid) );
    dumpTimes(i) = [sscanf(ln,'%f')]';
    i = i+1;
 end
fclose(fid);

%dumpTimes = dumpTimes(1:end-1); % sometimes the last dump time is missing

%---------------- Write file

fileName = strcat('gridSizePdf.dat');
fid = fopen(fileName, 'w');
fprintf(fid, '# times -->\n# ');
fprintf(fid, '%f  ', dumpTimes);
fprintf(fid, '\n# y_(m), Pdf...\n');
data = [xbins Pdx];
[ni nj] = size(data);
for i=1:ni
    fprintf(fid, '%-16.8e', data(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

%------------- plot output

ymin = 0.0;
ymax = max(max(Pdx))*1.1;

hgexport(gcf, 'gridSizePdf.pdf', hgexport('factorystyle'), 'Format', 'pdf');

hFig = figure('visible', 'off');
clf;

[T,X]=meshgrid(dumpTimes(2:end),xbins);
contourf(X,T,Pdx(:,2:end),50);
%contourf(X,T,Pdx(:,1:end),50);
shading flat;
title('Cell Size PDF', 'FontSize', 16);
xlabel('log_{10}({\Delta}x/(m))', 'FontSize', 16);
ylabel('Time (s)', 'FontSize', 16);
set(gca,'FontSize',16);
h = colorbar('East', 'FontSize', 16);
set(h,'YColor',[1 1 1]);
L=line([dxmin,dxmin], [dumpTimes(2),dumpTimes(end)]);
set(L,'Color',[1 1 1]);

hgexport(gcf, 'gridSizePdf.pdf', hgexport('factorystyle'), 'Format', 'pdf');

exit;

