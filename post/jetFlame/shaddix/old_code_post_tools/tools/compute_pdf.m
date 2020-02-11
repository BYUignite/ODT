
function [xbins P] = compute_pdf(x, nbins, xmin, xmax)

    if(nargin==1)
        nbins = 100;
        xmin = min(x);
        xmax = max(x);
    end
    if(nargin==2)
        xmin = min(x);
        xmax = max(x);
    end
    if(nargin==3)
        error('wrong number of arguments');
    end

    dxbin = (xmax-xmin)/nbins;
    xbins = linspace(xmin+dxbin/2, xmax-dxbin/2, nbins)';

    ibin = floor((x-xmin)/dxbin) + 1;
    ibin(find(ibin<1)) = 1;
    ibin(find(ibin>nbins)) = nbins;

    for i=1:nbins
        P(i) = length(find(ibin==i));
    end

    P = P/length(x)/dxbin;

end

