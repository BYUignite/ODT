
function w = fwhm(x,y)

   m = max(y);
   hm = m/2;
   im = find(y==m);

   ir = find(y>0.1*m & y<hm);
   il = min(ir);
   ih = max(ir);

   %x1 = interp1(y(1:im),x(1:im),hm,'linear', 'extrap');
   %x2 = interp1(y(im+1:end),x(im+1:end),hm,'linear', 'extrap');
   x1 = interp1(y(il:im),x(il:im),hm,'linear', 'extrap');
   x2 = interp1(y(im+1:ih),x(im+1:ih),hm,'linear', 'extrap');

   w = x2-x1;

end

