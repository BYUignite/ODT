% run as standalone (uncomment as noted), or use driver.sh
%
% Computes the PDF of accepted eddies
% Read in a file called out that has sizes of each accepted eddy as the columnPos column.
% Change nbins, and columnPos
% Plots the pdf with the model pdf and the model pdf with a new beta
%
% CHANGE SAMPLED PDF parameters below.


nbins = 60;
columnPos = 1;

%---------------- compute eddy pdf

oo = load('eddySizes.dat');

dat = oo(:,columnPos);
[npts id] = size(dat);

Lmin = 0.05*min(dat);
Lmax = 2.3*max(dat);
Lmin = 0.001;
Lmax = 10;
Lmin = min(dat);
Lmax = max(dat);

binf = logspace(log10(Lmin), log10(Lmax), nbins+1);
pdf = zeros(nbins,1);

for i=1:nbins

    jpos = find( dat >= binf(i) & dat < binf(i+1) );
    [ni id] = size(jpos);
    pdf(i) = pdf(i) + ni;

end

pdf = pdf / npts;
ds = zeros(nbins,1);
for i=1:nbins
    ds(i) = binf(i+1) - binf(i);
    pdf(i) = pdf(i) / ds(i);
end
ds = ds';

i=1:nbins;
bins = (binf(i+1) + binf(i))/2;

%------------- compare to sampled pdf

domainLength = 0.015232;
%domainLength = 0.017408;
%domainLength = 0.019278;
lp   = 0.015 * domainLength;
lmax = 1.0   * domainLength;
lmin = 0.004 * domainLength;
Cmax = exp(-2*lp/lmax);
Cmin = exp(-2*lp/lmin);
C    = 2*lp/(Cmax-Cmin);
f    = C./bins.^2 .*exp(-2*lp./bins);
f = f';

%------------- save results

data = [bins' pdf f];
save 'eddyPDF.dat' data;

exit



