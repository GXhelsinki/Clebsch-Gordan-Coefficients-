
% Code purpose:
% Compare the computational accuracy of Clebsch-Gordan coefficients;
% i.e., m2 range from m2_a to m2_b
%      m2_a <= m2 <= m2_b ;
%
% References:

% [1] Improving the recursion computation of Clebsch-Gordan coefficients for
% light scattering simulations. G. Xu. submitted to Journal of Quantitative
% Spectroscopy and Radiative Transfer Jun. 2020.

% [2] K. Schulten, R. Gordon, Recursive evaluation of 3j and 6j coefficients,
% Computer Physics Communications.

%Copyright@2020 Guanglang Xu.

%email : guanglang.xu@helsinki.fi;


clear; close all; clc;

j123 = [280, 220, 189];
m123 = [90, 135, 169 ];

n1 = j123(1);
n2 = j123(2);
n3=j123(3);
m1 = m123(1);
m2 = m123(2);

% scaling factor for obtaining 3j sum
j2tj3=(2*n3+1)/(2*n1+1);
% range of m2;

m2a = -min(n2, n3+m1);
m2b = min(n2, n3-m1);

Nmax=max(n1+n2, n2);
Nmax=max(Nmax,n3)+20000000;


tik1 = cputime;
[cgn1,m2_min,m2_max] = CG_coefficients_m2_range(n1, n2, n3, m2a, m2b, m1, Nmax);
tod1 = (cputime-tik1)*10^10;
Sos1 = sum(cgn1.^2, 'all');
R1= abs(1-Sos1/j2tj3);
%


tik2 = cputime;
[cgn2, sne1, kne1] = CG_coefficients_m2_range_log(n1,...
    n2, n3, m2a, m2b, m1, Nmax);
tod2 = (cputime-tik2)*10^10;
Sos2= sum(cgn2.^2, 'all');
R2= abs(1-Sos2/j2tj3);

acg = [cgn1 cgn2]; maxy=max(acg); miny=min(acg);
rt=1.2;
% ploting

mt2 =max(m2a, m2b);
mt1= min(m2a, m2b);
Nm = mt2-mt1+1;

close all;
plot(mt1:mt2, cgn2, '.r');
hold on;
plot(mt1:mt2, cgn1, '.k');
%hold on;
xlim([mt1*rt mt2*rt]);
ylim([miny*rt maxy*rt]);
legend('sign-exponent', 'three-term linear', 'Interpreter','latex');
legend('location', 'best');
legend boxoff;
ylabel('CG coefficients', 'Interpreter','latex' );
xlabel('magnetic quantum number $m_{2}$','Interpreter','latex');
xlo=mt1+(mt2-mt1)*0.12; ylo=miny+(maxy-miny)*0.06;

j1s=strcat('$j1=', num2str(n1),'$,',' ');
j2s=strcat('$j2=',num2str(n2),'$, ', ' ');
j3s=strcat('$j3=', num2str(n3),'$, ', ' ');
m1s=strcat('$m1=', num2str(m1),'$, ', ' ');
m2s=strcat('$m2_{min}=', num2str(mt1),'$,', ' ');
m2b=strcat('$m2_{max}=', num2str(mt2), '$.', ' ');

quantum_number_string = strcat(j1s, j2s, j3s, m1s, m2s, m2b);
text(xlo,ylo, quantum_number_string, ...
    'Interpreter','latex');

% m3=m1+m2;
% ft=(-1)^(n1-n2+m3); n=n3p;
% w3j = cgn1./sqrt(2*n+1).*ft;












