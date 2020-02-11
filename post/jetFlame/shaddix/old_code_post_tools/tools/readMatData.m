function A = readMatData(fname)

   fid = fopen(fname, 'r');

   i = 1;
   while(~feof(fid))

       ln = strtrim( fgetl(fid) );
       if(numel(ln)==0 || ln(1) == '#')
           continue;
       else
           A(i,:) = [sscanf(ln,'%f')]';
           i = i+1;
       end

   end

  fclose(fid);

end



