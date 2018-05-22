function write_profile(name,x,y,n,endian)

% write profile in file 

  if nargin<5, endian = 'n'; end
  
  prec = 'real*8';
  
  fid = fopen([name '.curve' ],'w',endian);
  fwrite(fid,x,prec);
  fwrite(fid,y,prec);
  fwrite(fid,n(:,1),prec);
  fwrite(fid,n(:,2),prec);
  fclose(fid);
