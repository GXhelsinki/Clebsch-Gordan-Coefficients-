function  [cgn, sne, kne] = CG_coefficients_m2_range_log(j1, j2, j3, m2_a, m2_b, m1, N_max)

% Code purpose:
% To computate a particular range of Clebsch-Gordan coefficients;
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
% sne -----------------% signs of the CG coefficients; 
% kne -----------------% exponents of the CG coefficients; 

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

% initializa outputs
cgn(1:N_cg)=zeros;
sne(1:N_cg)=ones;
kne(1:N_cg)=zeros;


if (j3>Jmax)||(j3<Jmin)
    fprintf('j3 is not within the range.');
    return;
end
if (abs(m1)>j1)
    fprintf('m1 is not within the range.');
    return;
end
if (j1>N_max)||(j2>N_max)||(j3>N_max)
    fprintf('The principle quantum number is too large.');
    return;
end


m2_min = -min(j2, j3+m1);   % smallest m2 for non-zero CGs;
m2_max = min(j2, j3-m1);    % largest m2 for non-zero CGs;


dntop = 10;   % number of critical terms

dn1= m2_max - m2_min;
Nc=dn1+1;

%upb = max(m2_max, m2_b);   % upper bound;
%lwb = min(m2_min, m2_a);   % lower bound;

%N_max_cg = upb-lwb+1;

Cx(Nc)=zeros;
Sn(Nc)=ones;
Kn(Nc)=zeros;

% obtain the sign , exponent, and value;
% j1, j2, j3, -m1, -m2_min

[Cg_min ,Sn_min, Kn_min] = CG_coefficients_n3_log(j1, j2, j3, -m1, -m2_min);


m3=-m1-m2_min;

%  sign
Sn(1)=Sn_min*(-1)^(j3+m3+2*j1);

% exponent;
Kn(1)=Kn_min+log(1/sqrt(2*j3+1));

%
% Cx(1) = Cg_min*(-1)^(j3+m3+2*j1)/sqrt(2*j3+1);
% value
if (Cg_min==0)
    Cx(1)= Sn(1)*exp(Kn(1));
else
    Cx(1) = Cg_min*(-1)^(j3+m3+2*j1)/sqrt(2*j3+1);
end


if (dn1 == 0) % only one value;
    
elseif (dn1==1) % two values ;
    
    m2 = m2_min;
    m3 = -m1-m2;
    
    dd2 =Df(m2, m3, j1, j2, j3);
    cc2 =Cf(m2+1, j2, j3, m3-1);
    
    %
    tp = -dd2/cc2*Sn(1);
    Sn(2)=sign(tp);
    
    Kn(2)=Kn(1)+log(abs(tp));
    Cx(2)=Sn(2)*exp(Kn(2));
    
    
    %Cx(2) = -dd2*Cx(1)/cc2;
    %
    
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
            
            tp = -dd2/cc2*Sn(k);
            Sn(k+1)=sign(tp);
            
            Kn(k+1)=Kn(k)+log(abs(tp));
            Cx(k+1)=Sn(k+1)*exp(Kn(k+1));
            
            %Cx(k+1)=-Cx(k)*dd2/cc2;
            %sc=sc+Cx(k+1)^2;
            
            k=k+1;
        else
            cc1=Cf(m2, j2, j3, m3);
            cc2= Cf(m2+1, j2, j3, m3-1);
            dd2=Df(m2, m3, j1, j2, j3);
            
            if(Sn(k)~=0)&&(Sn(k-1)~=0)
                
                Dm=Kn(k)-Kn(k-1);
                tp=-(dd2*Sn(k)+cc1*Sn(k-1)*exp(-Dm))/cc2;
                Sn(k+1)=sign(tp);
                
                Kn(k+1)=Kn(k)+log(abs(tp));
                
                Cx(k+1)=Sn(k+1)*exp(Kn(k+1));
                %Cx(k+1)=-(dd2*Cx(k)+Cx(k-1)*cc1)/cc2;
                %sc = sc + Cx(k+1)^2;
                k =k+1;
            else
                
                Cx(k+1)=-(dd2*Cx(k)+Cx(k-1)*cc1)/cc2;
                Kn(k+1)=log(abs(Cx(k+1)));
                Sn(k+1)=sign(Cx(k+1));
                k=k+1;
                
            end
            
        end
        
    end
elseif (dn1>dntop) % should use both upward and downward ;
    
    k=1;
    mid_dn = round(Nc/2)+1; % mid-range;
    im2b = m2_min+mid_dn;
    % upward reccurence;
    
    for m2 = m2_min:im2b
        m3 = -m1-m2;
        
        if (m2==m2_min)
            %k=k+1;
            dd2=Df(m2, m3, j1, j2, j3);
            cc2=Cf(m2+1, j2, j3, m3-1);
            
            tp = -dd2/cc2*Sn(k);
            Sn(k+1)=sign(tp);
            
            Kn(k+1)=Kn(k)+log(abs(tp));
            Cx(k+1)=Sn(k+1)*exp(Kn(k+1));
            
            %Cx(k+1)=-Cx(k)*dd2/cc2;
            %sc=sc+Cx(k+1)^2;
            
            k=k+1;
        else
            cc1=Cf(m2, j2, j3, m3);
            cc2= Cf(m2+1, j2, j3, m3-1);
            dd2=Df(m2, m3, j1, j2, j3);
            
            
            if(Sn(k)~=0)&&(Sn(k-1)~=0)  % avoid zeros;
                Dm=Kn(k)-Kn(k-1);
                tp=-(dd2*Sn(k)+cc1*Sn(k-1)*exp(-Dm))/cc2;
                Sn(k+1)=sign(tp);
                Kn(k+1)=Kn(k)+log(abs(tp));
                Cx(k+1)=Sn(k+1)*exp(Kn(k+1));
                k =k+1;
                
            else
                
                Cx(k+1)=-(dd2*Cx(k)+Cx(k-1)*cc1)/cc2;
                Kn(k+1)=log(abs(Cx(k+1)));
                Sn(k+1)=sign(Cx(k+1));
                k=k+1;
                
            end
            
        end
        
    end
    % maximum m2;
    [Cg_max, Sn_max, Kn_max]= CG_coefficients_n3_log(j1,j2, j3, -m1, -m2_max);
    % converted it to 3j-symbols;
    m3=-m1-m2_max;
    Sn(Nc)=Sn_max*(-1)^(j3+m3+2*j1);
    Kn(Nc) = Kn_max + log(1/sqrt(2*j3+1));
    
    if(Cg_max==0)
        Cx(Nc)=Sn(Nc)*exp(Kn(Nc));
    else
        Cx(Nc)=Cg_max*(-1)^(j3+m3+2*j1)/sqrt(2*j3+1);
    end
    
    
    k=Nc;
    
    for m2=m2_max:-1:im2b+1
        
        m3 = -m1-m2;
        
        if (k==Nc)
            
            cc2=Cf(m2, j2, j3, m3);
            dd2=Df(m2, m3, j1, j2, j3);
            
            tp=-dd2*Sn(k)/cc2;
            Sn(k-1)=sign(tp);
            Kn(k-1)=Kn(k)+log(abs(tp));
            Cx(k-1)=Sn(k-1)*exp(Kn(k-1));
            %Cx(k-1)=-dd2*Cx(k)/cc2;
            %sc2 = sc2 + Cx2(k+1)^2;
            k = k-1;
            
        else
            
            cc2=Cf(m2+1, j2, j3, m3-1);
            cc1=Cf(m2, j2, j3, m3);
            dd2 = Df(m2, m3, j1, j2, j3);
            
            if(Sn(k+1)~=0)&&(Sn(k)~=0)  % advoid zeros;
                Dm=Kn(k)-Kn(k+1);
                tp=-(cc2*Sn(k+1)*exp(-Dm)+dd2*Sn(k))/cc1;
                Sn(k-1)=sign(tp);
                Kn(k-1)=Kn(k)+log(abs(tp));
                Cx(k-1)=Sn(k-1)*exp(Kn(k-1));
                k = k-1;
            else
                
                Cx(k-1) = -(cc2*Cx(k+1)+dd2*Cx(k))/cc1;
                Kn(k-1)=log(abs(Cx(k-1)));
                Sn(k-1)=sign(Cx(k-1));
                k=k-1;
                
            end
            
        end
        
    end
    
    
end


for ik = m2_a:m2_b
    if (ik<m2_min)
        cgn(ik-m2_a+1)=0.0;
    elseif (ik>m2_max)
        cgn(ik-m2_a+1)=0.0;
    else
        
        %cgn(ik-m2_a+1) = Cx(ik-m2_min+1);
        
        m3=m1+ik;
        sne(ik-m2_a+1) = Sn(ik-m2_min+1)*(-1)^(j1-j2+m3);
        kne(ik-m2_a+1) = Kn(ik-m2_min+1) + log(sqrt(2*j3+1));
        cgn(ik-m2_a+1) = sne(ik-m2_a+1)*exp(kne(ik-m2_a+1));
        
        %cgn(ik-m2_a+1) = Cx(ik-m2_min+1)*sqrt(2*j3+1)*(-1)^(j1-j2+m3);
    end
    
end


    function C_m2m1=Cf(m2, j2, j3, m3)
        C_m2m1 = (j2-m2+1)*(j2+m2)*(j3+m3+1)*(j3-m3);
        C_m2m1 = sqrt(C_m2m1);
    end

    function D_m2=Df(m2, m3, j1, j2, j3)
        D_m2 = j2*(j2+1)+j3*(j3+1)-j1*(j1+1)+2.0*m2*m3;
    end


end