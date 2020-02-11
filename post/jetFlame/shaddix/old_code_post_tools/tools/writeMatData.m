function dummy = writeMatData(fname, A)

   fid = fopen(fname, 'w');

   [ni nj] = size(A);

   for i=1:ni
       fprintf(fid, '%-16.8e', A(i,:));
       fprintf(fid, '\n');
   end

   fclose(fid);

end
