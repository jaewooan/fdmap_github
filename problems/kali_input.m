% script illustrating how to enter spatially variable friction parameters,
% prestress, and state for FDMAP simulations with rate-and-state friction

ny = 301; y = linspace(-30,0,ny)'; % coordinates along fault, column vector

a = 0.01+0.005*sin(2*pi*y/10);
b = 0.016+0.004*exp(-0.5*((y+15)/3).^2);
V0 = 1e-6+0*y;
f0 = 0.6+0*y;
L = 0.05+0*y;
fw = 0*y; % required, even if not used
Vw = 1e10+0*y; % required, even if not used

% initial state
Psi = 0.6+0.1*exp(y/5);

% prestress (in addition to resolved stresses in body)
S0 = cos(2*pi*y/10);
N0 = sin(2*pi*y/4);

% write files
filename = 'kali.friction'; % any name is ok
fdmap_write_ratestate_friction(filename,a,b,V0,f0,L,fw,Vw)
filename = 'kali.state'; % any name is ok
fdmap_write_state(filename,Psi)
filename = 'kali.prestress'; % any name is ok
fdmap_write_fault_prestress(filename,S0,N0);