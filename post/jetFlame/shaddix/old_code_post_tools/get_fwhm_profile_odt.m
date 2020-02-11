% Read list of files and compute the fwhm for each one

clc;
clear;


myfile = 'means_6.dat';                 % holds mixture fraction. if mixf is not in column 6, then change here

f1 = fopen(myfile);

i=1;
while(~feof(f1))
    ln = fgetl(f1);
    A(i,:) = [sscanf(ln,'%f')]';
    i = i+1;
end
fclose(f1);

[npts, ntim] = size(A);
ntim=ntim-1;

width=zeros(ntim,1);
peak=zeros(ntim,1);

for i=1:ntim
    i
    width(i) = fwhm(A(:,1), A(:,i+1));
    peak(i)  = A(floor(npts/2),i+1);
end

times = zeros(ntim,1);

f1 = fopen('../input/dumpTimes.inp');
ln = fgetl(f1);
i=1;
while(~feof(f1))
    ln = fgetl(f1);
    times(i) = [sscanf(ln,'%f')]';
    i = i+1;
end

data = [times(1:end) width peak];

save -ascii fwhm_odt.dat data;

exit;


