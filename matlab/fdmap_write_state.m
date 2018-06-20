function fdmap_write_ratestate_friction(name,Psi)
    
    % write file

    prec = 'real*8'; endian = 'n';

    fid = fopen(name,'w',endian);
    fwrite(fid,Psi,prec);
    fclose(fid);
