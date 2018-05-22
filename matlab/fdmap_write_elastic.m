function fdmap_write_elastic(name,rho,cs,cp,endian)

% write elastic properties in file 

  if nargin<5, endian = 'n'; end
      
  prec = 'real*8';
  
  [fid,m] = fopen(name,'w',endian);
  if fid==-1,disp(m),return,end
  fwrite(fid,rho,prec);
  fwrite(fid,cs ,prec);
  fwrite(fid,cp ,prec);
  fclose(fid);
  

