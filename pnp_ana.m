function [C_pls, C_min, Ey, V] = pnp_ana(c, sig, dy, n_y)


%% constants

% k_B = 1.38e-23; % J/K
k_B = 1.38e-16; % erg/K; 1 erg = 10^-7 J
eps_0 = 8.85e-12; % F/m
% c_e = 1.602e-19; % C in SI
c_e = 4.80320425e-10; % C in cgs

%% exp parameters

T = 300;
eps = 80;
% sigma = -c_e / 250e-16; % example 1
% sigma = -c_e / 25e-16; % example 2
sigma = -sig * c_e / 10000000e-16; % custom
% n0 = 1; %M % example 1
% n0 = 0.1; %M % example 2
n0 = c; %M % custom


%% Debye length

l_B = c_e^2 / (eps*k_B*T)
l_D = sqrt( 8*pi*l_B* n0*6e23/1000 )^-1 % cm

%% counter-ion only

% l_B = c_e^2 / (4*pi*eps_0*eps*k_B*T)

b = eps*k_B*T / (2*pi*c_e*abs(sigma))

% b = eps*eps_0*T;% Gouy-Champman length
% v_ci = 

z1 = 0:(100*dy):( (n_y-1)*100*dy );
z2 = z1;


%% salt



gamma = -b/l_D + sqrt( (b/l_D)^2 + 1 );

v_salt = -2*k_B*T/c_e * log( (1+gamma*exp(-z1./l_D)) ./ (1-gamma*exp(-z1./l_D)) );

n_pls = n0 * ( (1+gamma*exp(-z2./l_D)) ./ (1-gamma*exp(-z2./l_D)) ).^2;
n_min = n0 * ( (1-gamma*exp(-z2./l_D)) ./ (1+gamma*exp(-z2./l_D)) ).^2;

%% generate C and P profiles for numerics

Ey = (-1/dy) * 299.792458 * diff(v_salt);


%% ouput

V = 299.792458*v_salt;
C_pls = n_pls;
C_min = n_min;

end







