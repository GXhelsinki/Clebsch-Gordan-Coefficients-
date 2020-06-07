function  [cgn, m2_min, m2_max] = CG_coefficients_m2_scaling(j1, j2, j3, m2_a, m2_b, m1, N_max)

% Code purpose:
% To computate a particular range of Clebsch-Gordan coefficients using the
% three-term linear recursion method by Schulten, Gordon wiht scaling and
% normailization;
% i.e., m2 range from m2_a to m2_b
%      m2_a <= m2 <= m2_b ;
%
% Input and Output
% ( j1, j2, j3, m2a , m2b, m1, Nmax)--->[cgn, sne, kne]

% Inputs:
% j1, j2, j3; -------------% principle quantum numbers
% m1  -----------------% magnetic quantum numbers for j1
% m2a ---------------- % minimum of m2 :
% m2b ---------------- % maximum of m2 ;
% Nmax ----------------% the upper limit of principle quantum numbers.

%Outputs:
% cgn -----------------% CGn coefficients ;
% m2_min -----------------% possible miminum m2;
% m2_nax -----------------% maximum m2;

% References:

% [1] M. I. Mishchenko, L. D. Travis, A. A. Lacis, Scattering, absorption,
% and emission of light by small particles (Appendix. D), Cambridge university press, 2002

% [2] Improving the recursion computation of Clebsch-Gordan coefficients for
% light scattering simulations. G. Xu. submitted to Journal of Quantitative
% Spectroscopy and Radiative Transfer Jun. 2020.

% [3] K. Schulten, R. Gordon, Recursive evaluation of 3j and 6j coefficients,
% Computer Physics Communications.

%Copyright@2020 Guanglang Xu.

%email : guanglang.xu@helsinki.fi;

% exlucdes the zero CG-coefficients;

Jmax = j1+j2; Jmin = abs(j1-j2);

if (m2_a>m2_b)
    m2_t =m2_a;
    m2_a =m2_b;
    m2_b =m2_t;
end


N_cg = m2_b-m2_a+1 ; % total number of CGs requested;

cgn(1:N_cg)=zeros;

if (j3>Jmax)||(j3<Jmin)
    return;
end
if (abs(m1)>j1)
    return;
end
if (j1>N_max)||(j2>N_max)||(j3>N_max)
    return;
end

%NP_L1= 2*N_max +1;         % possible number of non-zero coefficients;

m2_min = -min(j2, j3+m1);   % smallest m2 for non-zero CGs;
m2_max = min(j2, j3-m1);    % largest m2 for non-zero CGs;

% kp = max(0, m2_min-m2_a);
% if (kp<0), m2_min<m2_a; m2_a>m2_min;

% kp2 = max(0, m2_b-m2_max);

dntop = 10;

% im2 = m2_min;
% cm2 = m2_max;

dn1= m2_max - m2_min;

Nc = dn1+1;

upb = max(m2_max, m2_b);   % upper bound;
lwb = min(m2_min, m2_a);   % lower bound;

%N_max_cg = upb-lwb+1;

Cx(Nc)=zeros;
Cx2=Cx;

Cx(1) = 1.0; % starting with 1.0;
%cmax=Cx(1);
sc = Cx(1)*Cx(1);


if (dn1 == 0) % only one value;
elseif (dn1==1) % two values ;
    
    m2 = m2_min;
    m3 = -m1-m2;
    
    dd2 =Df(m2, m3, j1, j2, j3);
    cc2 =Cf(m2+1, j2, j3, m3-1);
    
    Cx(2) = -dd2*Cx(1)/cc2;
    
    sc = sc + Cx(2)^2;
    
    %sc = 1.0/(sqrt(2*j1+1))*1.0/sqrt(Cx(1)^2 +Cx(2)^2)*sign(Cx(2))*...
    %    (-1)^(j2-j3-m1);
    
    % Cx(1) = Cx(1)*sc;
    % Cx(2) = Cx(2)*sc;
    
elseif (dn1<dntop) % upward reccurence ;
    
    k=1;
    %m2 = m2_min;
    for m2 = m2_min:m2_max
        %k=k+1;
        m3 = -m1 - m2;
        if (m2==m2_min)
            %k=k+1;
            dd2=Df(m2, m3, j1, j2, j3);
            cc2=Cf(m2+1, j2, j3, m3-1);
            Cx(k+1)=-Cx(k)*dd2/cc2;
            sc=sc+Cx(k+1)^2;
            
            k=k+1;
        else
            cc1=Cf(m2, j2, j3, m3);
            cc2= Cf(m2+1, j2, j3, m3-1);
            dd2=Df(m2, m3, j1, j2, j3);
            Cx(k+1)=-(dd2*Cx(k)+Cx(k-1)*cc1)/cc2;
            sc = sc + Cx(k+1)^2;
            k =k+1;
            
        end
        
    end
elseif (dn1>dntop) % should use both upward and downward ;
    
    k=1;
    
    mid_dn = round(dntop/2)+1; % mid-range;
    
    im2b = m2_min+mid_dn;
    
    
    
    % upward reccurence;
    
    for m2 = m2_min:im2b+1
        m3 = -m1-m2;
        
        if (m2==m2_min)
            %k=k+1;
            dd2=Df(m2, m3, j1, j2, j3);
            cc2=Cf(m2+1, j2, j3, m3-1);
            Cx(k+1)=-Cx(k)*dd2/cc2;
            sc=sc+Cx(k+1)^2;
            
            k=k+1;
        else
            cc1=Cf(m2, j2, j3, m3);
            cc2= Cf(m2+1, j2, j3, m3-1);
            dd2=Df(m2, m3, j1, j2, j3);
            Cx(k+1)=-(dd2*Cx(k)+Cx(k-1)*cc1)/cc2;
            sc = sc + Cx(k+1)^2;
            k =k+1;
            
        end
        
    end
    
    N_upward= im2b+1-m2_min+1;  % number of upward;
    %N_downward = m2_max-(im2b-1)+1; % number of downward;
    
    % downward;
    k = Nc;
    Cx2(Nc)=1.0;
    sc2= Cx2(Nc)^2;
    
    for m2=m2_max:-1:im2b-1
        
        m3 = -m1-m2;
        
        if (m2==m2_max)
            cc2=Cf(m2, j2, j3, m3);
            dd2=Df(m2, m3, j1, j2, j3);
            
            Cx2(k-1)=-dd2*Cx2(k)/cc2;
            sc2 = sc2 + Cx2(k-1)^2;
            k = k-1;
            
        else
            
            cc2=Cf(m2+1, j2, j3, m3-1);
            cc1=Cf(m2, j2, j3, m3);
            dd2 = Df(m2, m3, j1, j2, j3);
            Cx2(k-1) = -(cc2*Cx2(k+1)+dd2*Cx2(k))/cc1;
            sc2=sc2 + Cx2(k-1)^2;
            k = k-1;
            
        end
        
    end
    
    % index of the same number;
    m_int_up = N_upward-1;
    m_int_down = m_int_up;
    
    up = Cx(m_int_up+1)*Cx2(m_int_down+1)+Cx(m_int_up)*Cx2(m_int_down)...
        +Cx(m_int_up-1)*Cx2(m_int_down-1);
    
    dn = Cx(m_int_up+1)^2 + Cx(m_int_up)^2 + Cx(m_int_up-1)^2;
    
    lamda = up/dn;
    
    %lamda = Cx2(m_int_down)/Cx(m_int_up);
    
    Cx=Cx*lamda;
    
    sc=sc*lamda^2;
    
    sc=sc+sc2-Cx(m_int_up+1)^2-Cx2(m_int_down-1)^2-Cx(m_int_up)^2;
    
    %mtt = m_int_up+m_int_down;
    
    Cx(m_int_up+1: end) = Cx2(m_int_up+1: end);
    
    %Cx=[Cx(1:m_int_up) flip(Cx2(1:m_int_down))];
    
    
end

% now we have the non-scaled 3j symbols Cx and sc, the summation of squares;

% kl = length(Cx)

sc = 1/(2*j1+1)/sc;
sc = sqrt(sc);
%im2_max = dn1+1;

sc = sign(Cx(Nc))*(-1)^(j2-j3-m1)*sc ;

for ik = m2_a:m2_b
    if (ik<m2_min)
        cgn(ik-m2_a+1)=0.0;
    elseif (ik>m2_max)
        cgn(ik-m2_a+1)=0.0;
    else
        cgn(ik-m2_a+1) = Cx(ik-m2_min+1) * sc;
        cgn(ik-m2_a+1) = cgn(ik-m2_a+1)*sqrt(2*j3+1)*(-1)^(j1-j2+m1+ik);
        
    end
    
end


    function C_m2=Cf(m2, j2, j3, m3)
        C_m2 = (j2-m2+1)*(j2+m2)*(j3+m3+1)*(j3-m3);
        C_m2 = sqrt(C_m2);
    end
    function D_m2=Df(m2, m3, j1, j2, j3)
        D_m2 = j2*(j2+1)+j3*(j3+1)-j1*(j1+1)+2.0*m2*m3;
    end


end