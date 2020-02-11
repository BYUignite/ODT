
path(path,'../tools');

data = readMatData('eddyData.dat');
Leddy = data(:,2);
logLeddy = log10(Leddy);

[xbins_log Peddy_log] = compute_pdf(logLeddy, 80);

command = '../tools/getInputFileParameter.py ../../input/odtParam.inp domainLength';
[status domainLength] = system(command);
domainLength = str2num(domainLength);
command = '../tools/getInputFileParameter.py ../../input/odtParam.inp Lmin';
[status Lmin] = system(command);
Lmin = str2num(Lmin);
Lmin = Lmin * domainLength;

%------------- plot output

ymin = min(Peddy_log);
ymax = max(Peddy_log)*1.1;
xmin = min(xbins_log); xmin = min([xmin,log10(Lmin)])-0.5;
xmax = max(xbins_log); xmax = max([xmax,log10(Lmin)])+0.5;

hFig = figure('visible', 'off');
clf;
plot(xbins_log, Peddy_log, [log10(Lmin),log10(Lmin)],[ymin,ymax]);
xlabel('log_{10}(L_{eddy})', 'FontSize', 16);
ylabel('P(log_{10}(L_{eddy}))', 'FontSize', 16);
ylim([ymin,ymax]);
xlim([xmin,xmax]);
set(gca,'FontSize',16);

hgexport(gcf, 'eddySizePdf.pdf', hgexport('factorystyle'), 'Format', 'pdf');

exit;

