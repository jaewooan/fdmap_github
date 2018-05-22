function fdmap_write_acoustic_gravity(name,drho0dy,endian)

% write acoustic gravity wave properties in file 

  if nargin<3, endian = 'n'; end
      
  prec = 'real*8';
  
  [fid,m] = fopen(name,'w',endian);
  if fid==-1,disp(m),return,end
  fwrite(fid,drho0dy,prec);
  fclose(fid);
  
