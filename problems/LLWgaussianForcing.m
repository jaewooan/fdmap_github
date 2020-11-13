% script to generate forcing file for LLW with Gaussian,
% for comparison to (partially) hard-coded forcing

% forcing filename (can be anything)
filename = 'LLWgaussian.forcing';

% forcing function

gauss_amp = 1;
gauss_sigma = 1.25;
gauss_sigmat = 0.125;

% anonymous function to faciliate evaluation
gauss = @(x,y,t) gauss_amp*exp(-0.5*(x/gauss_sigma).^2 - 0.5*(y/gauss_sigma).^2) .* exp(-0.5*((t-4*gauss_sigmat)/gauss_sigmat).^2) / (gauss_sigmat*sqrt(2*pi));

endian = 'n'; % native endian
prec = 'real*8'; % double precision

% open file
[fid,m] = fopen(filename,'w',endian);
if fid==-1,disp(m),return,end

% evaluate on grid and write to disk

% spatial grid
X = [-50:50]; nx = length(X);
Y = [-50:50]; ny = length(Y);
[x,y] = ndgrid(X,Y);

% arrays containing fields to be written
vz = zeros(nx,ny);
sxz = zeros(nx,ny);
syz = zeros(nx,ny);

% time
dt = 1.25; t = [0:100]*dt; nt = length(t); % note nt is inclusive of t=0

for n=1:nt % loop over time steps
  vz = gauss(x,y,t(n)); % evaluates forcing function
  pcolored(X,Y,vz) % check by plotting (optional)
  drawnow
  % write to disk
  fwrite(fid,vz,prec); % nonzero forcing
  fwrite(fid,sxz,prec); % must write these, too
  fwrite(fid,syz,prec); % must write these, too
end

% close file
fclose(fid);
