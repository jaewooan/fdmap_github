% script to load FDMAP simulation data and make plots

% add FDMAP matlab directory to your path

addpath ../matlab

% load basic simulation information
% first argument is problem name (prefix in name='...' in input file)
% second argument is directory in which files were written

A = init('gp229sw3ex1','../data');

% plot fields on fault

% coordinate vectors on fault
% only one index, which is grid point index along fault

X = loadfast(A,'I_x'); % x (km)
Y = loadfast(A,'I_y'); % y (km)

% fields, stored as data structures

D = loadfast(A,'I_D'); % slip (m)
V = loadfast(A,'I_V'); % slip velocity (m/s)
S = loadfast(A,'I_S'); % shear stress (MPa)
N = loadfast(A,'I_N'); % normal stress (MPa)

% loop over time steps, plot slip velocity
% can be done for other fields, too

figure(1)
for n=1:V.nt
    plot(X,V.V(:,n))
    xlabel('x (km)')
    ylabel('V (m/s)')
    title(['t = ' num2str(V.t(n)) ' s'])
    drawnow
end

% contours of slip every 10 time steps

figure(2)
plot(X,D.D(:,1:10:end),'b')
xlabel('x (km)')
ylabel('slip (m)')
title(['plotted every ' num2str(D.t(11)) ' s'])

% space-time plot of slip velocity
% can do this with other fields, too

figure(3)
pcolored(X,V.t,V.V)
xlabel('x (km)')
ylabel('t (s)')
title('slip velocity (m/s)')

% add line correponding to S-wave speed
cs = A.M{1}.cs; % S-wave speed (km/s)
hold on
plot([0 40],[0 40]/cs','w')
hold on

% load coordinate arrays in body
% these have two indices x(i,j) is at x coordinate value at
% grid point having indices (i,j)

x = loadfast(A,'B_x'); % x (km)
y = loadfast(A,'B_y'); % y (km)

% load particle velocity, stored as data structure
% vz.t = time vector
% vz.nt = number of time steps
% vz.vz = particle velocity, first two indices are space, third is time

vz = loadfast(A,'B_vz'); % v_z (m/s)

% loop over time steps and plot velocity

figure(4)
cmap % set blue-white-red colormap
for n=1:vz.nt
    pcolored(x,y,vz.vz(:,:,n))
    xlabel('x (km)')
    ylabel('y (km)')
    title('particle velocity v_z (m/s)')
    caxis([-5 5])
    drawnow
end

% you can do the same for the two stress components, loaded as

sxz = loadfast(A,'B_sxz'); % sigma_{xz} (MPa)
syz = loadfast(A,'B_syz'); % sigma_{yz} (MPa)
