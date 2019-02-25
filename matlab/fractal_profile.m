% self-similar fractal fault profile
% Eric M. Dunham
% original: November 13, 2012
% modified: February 25, 2019 (added comments)

% fault length L, grid spacing h=L/N
% periodicity is assumed for profile, so y(L/2)=y(-L/2)

% NOTE: fault will have ODD number N+1 grid points located at nodes
% x = [-L/2:h:L/2] (length L, grid spacing h=L/N)
% Fourier method will give N points with replica length L
% so an extra point, at same (x,y) as first point, is appended at end

write_files = false; % write files in FDMAP format
endian = 'n'; % endian of binary files

r = 1; % mesh refinement factor
N = 1000*r; % number of grid points minus one
h = 0.1/r; % grid spacing

alpha = 10^(-2); % amplitude-to-wavelength ratio of roughness
cutoff = true; Lmin = 20*h; % minimum wavelength (nominally 20*h)

L = N*h; % profile length
x = [-L/2:h:L/2]; % coordinate at nodes

shift_profile = true; % shift profile in y direction so that yend=0
x_shift = 0; % shift profile in x direction by this amount

% wavenumber
  
kN = pi/h; % Nyquist wavenumber
k = zeros(1,N);
k(1:N/2+1) = 2*[0:N/2]/N;
k(N:-1:N/2+2) = -k(2:N/2);
k = k*kN;

% white noise, unit normal distribution

seed = 17; % seed, set to specific value to get exactly same profile
s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);

y = randn([1,N]); % Gaussian white noise

% scale so PSD has unit amplitude

y = y*sqrt(N/L);

% FFT

Y = fft(y)*h;

% calculate PSD and check for unit amplitude

PSDy = abs(Y).^2/L;
disp(['PSD = ' num2str(mean(PSDy))])

% multiply Y by square-root of desired PSD

PSDy_exact = (2*pi)^3*alpha^2*k.^(-3);
Y = Y.*sqrt(PSDy_exact);

% remove k=0 component

Y(1) = 0;

% add short wavelength (high wavenumber) cutoff

if cutoff
  kmax = 2*pi/Lmin;
  I = find(abs(k)>kmax); 
  Y(I) = 0;
end

% inverse FFT, explointing conjugate symmetry of real-valued function

y = ifft(Y,'symmetric')/h;

% calculate fault slope

k(N/2+1) = 0; % no Nyquist for spectral differentiation
M = i*k.*Y;
m = ifft(M,'symmetric')/h;
  
% convert slope into unit normals

n = zeros(N,2);
n(:,1) = -m./sqrt(1+m.^2);
n(:,2) = 1./sqrt(1+m.^2);

% check alpha

alpha_check = sqrt(h*sum(y.^2)/L)/L;
ratio = alpha/alpha_check;
disp(['input alpha=',num2str(alpha),' calculated alpha=',num2str(alpha_check)])
disp(['ratio=',num2str(ratio)])

% only return non-negative portion of spectrum

k = k(1:N/2+1); k(N/2+1) = kN;
Y = Y(1:N/2+1);
M = M(1:N/2+1);
PSDy = abs(Y).^2/L;
PSDy_exact = PSDy_exact(1:N/2+1);

%return

% repeat first point, shift if needed

y = [y y(1)];
n = [n; n(1,:)];
x = x-x_shift;

xend = [-L/2 L/2]-x_shift; % endpoints of fault
yend = y(1); if shift_profile, y = y-yend; yend = 0; end

% plot profile (with vertical exaggeration)

clf,subplot(2,1,1)
plot(x,y,'b-',xend,yend,'ko')
disp('(not to scale)')

% rotate profile for dipping faults, given unit normal of nominal surface

nhat = [0 -1]; % unit normal of nominal surface
[n(:,1),n(:,2)] = rotate_nt2xy_vec(n(:,1),n(:,2),nhat);

% plot profile to scale

subplot(2,1,2)
plot(x,y,'b-',xend,yend,'ko')
hold on,quiver(x,y,n(:,1)',n(:,2)',0,'r'),hold off
axis image

if write_files % write files

  disp('USER CAN SET DIFFERENT FILE NAME'),return
  name = ['~/fdmap/profiles/alpha' sprintf('%2.0f',alpha) 'N' num2str(N+1) 'seed' num2str(seed)];
  fdmap_write_profile(name,x,y,n,endian)
  disp(['wrote profile ' name])

end

return

% power spectral density estimators

S = spectrum.mtm;
Hpsd = psd(S,y,'Fs',1/h,'SpectrumType','twosided');
loglog(k/(2*pi),PSDy,'r',Hpsd.Frequencies,Hpsd.Data,'b',k/(2*pi),PSDy_exact,'k'),xlim([1e-3 max(k/(2*pi))])

