
% Code purpose:
% Compute particular Clebsch-Gordan coefficient using
% the three-term linear method and sign-exponent recurrence method;;
%
% References:

% [1] Improving the recursion computation of Clebsch-Gordan coefficients for
% light scattering simulations. G. Xu. submitted to Journal of Quantitative
% Spectroscopy and Radiative Transfer Jun. 2020.

% [2] K. Schulten, R. Gordon, Recursive evaluation of 3j and 6j coefficients,
% Computer Physics Communications.

%Copyright@2020 Guanglang Xu.
%email : guanglang.xu@helsinki.fi;

%
clear; close all; clc;

j123 = [900, 620, 780];
m123 = [300, 350, 169 ];

n1 = j123(1);
n2 = j123(2);
n3=j123(3);
m1 = m123(1);
m2 = m123(2);

% scaling factor for obtaining 3j sum
j2tj3=(2*n3+1)/(2*n1+1);
% range of m2;
% specify the value of m2;
m2a=m2; m2b=m2;

%---------
Nmax=max(n1+n2, n2);
Nmax=max(Nmax,n3)+20000000;

tik1 = cputime;
[cgn1,m2_min,m2_max] = CG_coefficients_m2_range(n1, n2, n3, m2a, m2b, m1, Nmax);
tik1a=cputime;
tod1 = tik1a-tik1;

m3=m1+m2;
ft=(-1)^(n1-n2+m3); %n=n3;
w3j_linear = cgn1./sqrt(2*n3+1).*ft;
fprintf('%.5f\n',tod1);
fprintf('%.15f\n',w3j_linear);
fprintf('%.15f\n', cgn1);

tik2 = cputime;
[cgn2, sne1, kne1] = CG_coefficients_m2_range_log(n1,...
    n2, n3, m2a, m2b, m1, Nmax);
tik3=cputime;
tod2 = tik3-tik2;
w3j_sign_exp = cgn2./sqrt(2*n3+1).*ft;
fprintf('%.5f\n',tod2);
fprintf('%.15f\n',w3j_sign_exp);
fprintf('%.15f\n',cgn2);



















