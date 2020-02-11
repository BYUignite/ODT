%-------------------------------------------------------------------------------------
% DOL 12/28/09
% Script computes means, rms, conditional means, rms for all realizations in data folders at
% each of the dump times.
% Outputs data files like means_6.dat which is the mean for variable in column 6 of the data files.
% Also outputs chiStoic profiles.
%
% JCH note: weighting for conditional PDF needs to be multiplied by dx
%-------------------------------------------------------------------------------------

%------------ CHANGE THESE

%ntimes = 09;                            % # of times dumped for each realization **PASSED AS ARGUMENT** (set for standalone)
%ntimes = ntimes - 1;
%nrlz = 3;
Ldom   = 0.5;                            % domain size
nvar   = 38;                             % # of variables (columns in data file)
fstoic = 0.0637;
imixf  = 6;                              % column holding mixture fraction in data files
iposf  = 2;                              % column holding face position in data files
irho   = 3;                              % column holding density in data files
iuvel  = 11;
nx     = 501;                            % # of uniform grid points to interpolate to
nfbins = 60;                             % # of mixture fraction bins
ichi   = 7;                              % column holding chi in data files
Favre = 0;                               % set to true to do Favre averages

%------------

cmbl   = 1;                              % set to 0 to do the odtl's
X      = linspace(0,Ldom, nx); X = X';   % interpolated domain grid
df     = 1/nfbins;                       % mixf bin spacing
fbins  = linspace(df/2,1-df/2, nfbins);  % mixture fraction domain for conditional variables
fbins  = fbins';
varIndxList = 1:nvar;                    % set to the variable indicies you want
%varIndxList = [7];                      % set to the variable indicies you want

if(~exist('nrlz'))
    nrlz=128;
end
if(~exist('caseD'))
    error('missing caseD argument');
end
if(~exist('ntimes'))
    error('missing ntimes argument');
end

%-------------------------------------------------------------------------------------

means = zeros(ntimes, nvar, nx);         % processed data
mean2 = zeros(ntimes, nvar, nx);         % mean squares

rhouu = zeros(ntimes,1);
rhou  = zeros(ntimes,1);

cmeans = zeros(ntimes, nvar, nfbins);    % conditional means
cmean2 = zeros(ntimes, nvar, nfbins);    % conditional means
nBinHits = zeros(ntimes, nvar, nfbins);  % number of points in each mixf bin for cond means

%-------------------------------------------------------------------------------------

for itime = 1:ntimes                     % loop over dump times
    for ifile = 1:nrlz                   % loop over realizations

        disp(['Processing dump time #',num2str(itime),' for realization ',num2str(ifile)]);

        %if(ifile==113 | ifile==29 | ifile==47 | ifile==49 | ifile==77 | ifile==90)
        %    continue;
        %end

        if(cmbl == 1)
            command = strcat( 'ls ../data/', caseD, '/data_', num2str(ifile-1), ...
                              '/dmp_odtl_', num2str(itime),'.*' );
            %command = strcat( 'ls ../data/', caseD, '/data_', num2str(ifile-1), ...
            %                  '/dmp_cmbl_', num2str(itime),'_*' );
       else
            command = strcat( 'ls ../data/', caseD, '/data_', num2str(ifile-1), ...
                              '/dmp_odtl_', num2str(itime),'_*' );
        end

        [status myfile] = system(command);
        file = fopen(myfile(1:end-1));
        if(file==-1)
            error('ERROR could not open file ');
        end

        ln = fgetl(file);
        ln = fgetl(file);    % get the header lines
        ln = fgetl(file);
        if(cmbl==1) ln = fgetl(file); end

        clear A;                               % holds the contents of the file
        i = 1;
        while(~feof(file))                     % file A line by line in file
            ln = fgetl(file);
            A(i,:) = [sscanf(ln,'%f')]';
            i = i+1;
        end
        x = A(:,1);                            % current position array
        npts = length(x);

        % The weighting for statistics should be according to the volume
        % (or mass for Favre averaging).  Volume weighting does not need to
        % be done for the spatial averages that are interpolated onto a
        % uniform grid but should be done by density for Favre averaged.
        dvol = zeros(npts,1);
        wt   = zeros(npts,1);
        for ipt=1:npts-1
            dvol(ipt) = A(ipt+1,iposf) - A(ipt,iposf) ;
            if ( Favre )
                wt(ipt) = dvol(ipt) * A(ipt,irho);
            else
                wt(ipt) = dvol(ipt);
            end
        end

        % interpolate density to use in Favre averages on interpolated X grid
        % same as below in for ivar=1:nvar loop

        y = A(:,irho);
        rho = interp1(x,y,X,'spline','extrap');
        y = A(:,iuvel);
        uvel = interp1(x,y,X,'spline','extrap');
        rhouu(itime) = rhouu(itime) + sum( rho.*uvel.^2 );
        rhou(itime)  = rhou(itime)  + sum( rho.*uvel    );

        % for means tensors, indexing is means(itime,ivar,position)
        % for cmeans tensors, indexing is cmeans(itime,ivar,ibin)

        for ivar=1:nvar                        % add each variable to the "means"
            y = A(:,varIndxList(ivar));
            Y = interp1(x,y,X,'spline','extrap');
            Y2 = Y.*Y;

            if ( Favre )
                means(itime, ivar, :) = means(itime, ivar, :) + reshape(Y,1,1,nx) .* rho;
                mean2(itime, ivar, :) = mean2(itime, ivar, :) + reshape(Y2,1,1,nx) .* rho;
            else
                means(itime, ivar, :) = means(itime, ivar, :) + reshape(Y,1,1,nx);
                mean2(itime, ivar, :) = mean2(itime, ivar, :) + reshape(Y2,1,1,nx);
            end

            for ipt=1:npts
                ibin = floor(A(ipt,imixf)/df)+1;
                if(ibin <= 0) ibin = 1; end
                if(ibin >  nfbins) ibin = nfbins; end
                cmeans(itime, ivar, ibin) = cmeans(itime,ivar, ibin) + y(ipt) .* wt(ipt);
                cmean2(itime, ivar, ibin) = cmean2(itime,ivar, ibin) + y(ipt) .* y(ipt) .* wt(ipt);
                nBinHits(itime, ivar, ibin) = nBinHits(itime, ivar, ibin) + wt(ipt);
            end

        end

        fclose(file);

    end
end

uxt = rhouu./rhou;
save('uxt.dat', 'uxt', '-ascii');

means = means / nrlz;                          % normalize the means
mean2 = mean2 / nrlz;

sig2 = mean2 - means.*means;
sig = sqrt(abs(sig2));

cmeans = cmeans./nBinHits;
cmean2 = cmean2./nBinHits;
cmeans(find(nBinHits==0))=0;

csig2 = cmean2 - cmeans.*cmeans;
csig  = sqrt(abs(csig2));
csig(find(nBinHits==0))=0;


%-------------------------------------------------------------------------------------

chiStoic = cmeans(:,ichi,:);
ilo=max(find(fbins<=fstoic));
ihi=min(find(fbins>fstoic));
chi1 = cmeans(:,ichi,ilo);
chi2 = cmeans(:,ichi,ihi);
x1 = fbins(ilo);
x2 = fbins(ihi);
chiStoic = chi1+(chi2-chi1)/(x2-x1)*(fstoic-x1);

save('chiStoic.dat', 'chiStoic', '-ascii');

chiStoic = csig(:,ichi,:);
ilo=max(find(fbins<=fstoic));
ihi=min(find(fbins>fstoic));
chi1 = csig(:,ichi,ilo);
chi2 = csig(:,ichi,ihi);
x1 = fbins(ilo);
x2 = fbins(ihi);
chiStoic = chi1+(chi2-chi1)/(x2-x1)*(fstoic-x1);

save('chiStoicSig.dat', 'chiStoic', '-ascii');

%----------- save each variable in a data file:
%----------- first column is X, others are the var at the various times

for i=1:nvar
    var = means(:,i,:);
    var = reshape(var,ntimes, nx);
    var = [X var'];
    filename = strcat('means_', num2str(i), '.dat');
    save(filename, 'var', '-ascii');
end

for i=1:nvar
    var = sig(:,i,:);
    var = reshape(var,ntimes, nx);
    var = [X var'];
    filename = strcat('sig_', num2str(i), '.dat');
    save(filename, 'var', '-ascii');
end

for i=1:nvar
    var = cmeans(:,i,:);
    var = reshape(var,ntimes, nfbins);
    var = [fbins var'];
    filename = strcat('cmean_', num2str(i), '.dat');
    save(filename, 'var', '-ascii');
end

for i=1:nvar
    var = csig(:,i,:);
    var = reshape(var,ntimes, nfbins);
    var = [fbins var'];
    filename = strcat('csig_', num2str(i), '.dat');
    save(filename, 'var', '-ascii');
end

exit;
