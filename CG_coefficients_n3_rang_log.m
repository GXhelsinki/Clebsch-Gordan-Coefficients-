function  [cg, snt, knt] = CG_coefficients_n3_rang_log(n1, n2,n3a, n3b, m1, m2)

% Code Purpose:
% Compute a range of Clebsch-Gordan coefficient using
% the sign-exponent recurrence method;;

% Inputs:
% n1, n2 -------------% principle quantum numbers 
% n3a ----------------% requseted minimum n3;
% n3b ----------------% requested maximm  n3; 
% m1  ----------------% magnetic quantum numbers for j1
% m2 ---------------- % magnetic quantum numbers for j2

%Outputs:
% cgn -----------------% CGn coefficients ; 
% sne -----------------% signs of the CG coefficients; 
% kne -----------------% exponents of the CG coefficients; 

%
% References:

% [1] Improving the recursion computation of Clebsch-Gordan coefficients for
% light scattering simulations. G. Xu. submitted to Journal of Quantitative
% Spectroscopy and Radiative Transfer Jun. 2020.

% [2] K. Schulten, R. Gordon, Recursive evaluation of 3j and 6j coefficients,
% Computer Physics Communications.

%Copyright@2020 Guanglang Xu.
%email : guanglang.xu@helsinki.fi;


m3 = m2 + m1; % m' ;
N3min = max(abs(n1-n2), abs(m3));
N3min = max(N3min, 1);
N3max = n1+n2;

% exchange the order;

if (n3a>n3b)
    n3t =n3a;
    n3a =n3b;
    n3b =n3t;
end

% initialization;
No3=n3b-n3a+1;
cg=zeros(No3,1); 
snt=ones(No3,1);
knt=zeros(No3,1);

if (n3a>N3max)||(n3b<N3min)
    %cg=0.0;
    fprintf('j3 is not within the range.');
    return;
end

if (n1<0)||(n2<0)
    %cg=0.0;
    fprintf('n1 and n2 shoud be positive number.');
    return;
end
if (abs(m2)>n2) || (abs(m1)>n1)
    fprintf('m1 and m2 are not within the range.');
    %cg = 0.0;
    return;
end

if(N3max==0)
    fprintf('n3 max is zero.');
    %cg=0.0; 
    return;
end


% possible non-zeros;
nC = N3max-N3min+1;
Nmid = round((N3min + N3max)/2.);
CG(nC) = zeros; % CG-coefficients;
Sn(nC) = ones;  % sign of coefficients;
Kn(nC) = zeros; % exponents of the coefficients;

% masge1 =strcat('possible non-zero number is ', num2str(nC)); 
% masge2 =strcat('requested number of coefficients is ', num2str(No3)); 
% 
% fprintf(masge1);
% fprintf(masge2);


% kp1=n3p-N3min;
% kp2=N3max-n3p;

% np = Nmin;
if (( abs(n1-n2)>=abs(m3))&& (n1>=n2))
    u1 = logfact(n1+m1);
    u2 = logfact(n1-m1);
    u3 = logfact(2*n2);
    u4 = logfact(2*n1 - 2*n2 +1);
    d1 = logfact(2*n1+1);
    d2 = logfact(n2+m3-m1);
    d3 = logfact(n2-m3+m1);
    d4 = logfact(n1-n2+m3);
    d5 = logfact(n1-n2-m3);
    fact = (-1)^(n2+m3+m1);
    CG(1) = fact * exp((u1+u2+u3+u4-d1-d2-d3-d4-d5)/2);
    %ktt=u1+u2+u3+u4-d1-d2-d3-d4-d5
    %b=3
    Sn(1)=fact;
    Kn(1) = (u1+u2+u3+u4-d1-d2-d3-d4-d5)/2;
    CG(1)=Sn(1)*exp(Kn(1));
    
    
elseif  (( abs(n1-n2)>=abs(m3))&& (n1<n2))
    % mp == mp
    % mp-m == m ;
    % n1 == n;
    u1 = logfact(n2+m3-m1);
    u2 = logfact(n2-m3+m1);
    u3 = logfact(2*n1);
    u4 = logfact(2*n2 - 2*n1 +1);
    d1 = logfact(2*n2+1);
    d2 = logfact(n1+m1);
    d3 = logfact(n1-m1);
    d4 = logfact(n2-n1+m3);
    d5 = logfact(n2-n1-m3);
    fact = (-1)^(n1+m3+m3-m1);
    %CG(1) = fact * exp((u1+u2+u3+u4-d1-d2-d3-d4-d5)/2);
    
    Sn(1)=fact; % sign;
    Kn(1) = (u1+u2+u3+u4-d1-d2-d3-d4-d5)/2; % exponents;
    CG(1)=Sn(1)*exp(Kn(1));
    
    %ktt=u1+u2+u3+u4-d1-d2-d3-d4-d5
    %b=4
elseif (( abs(n1-n2)<abs(m3))&& (m3>=0))
    
    u1 = logfact(2*m3+1);
    u2 = logfact(n2+n1-m3);
    u3 = logfact(n1+m1);
    u4 = logfact(n2+m3-m1);
    d1 = logfact(n1+n2+m3+1);
    d2 = logfact(n1-n2+m3);
    d3 = logfact(n2-n1+m3);
    d4 = logfact(n1-m1);
    d5 = logfact(n2-m3+m1);
    fact = (-1)^(n1+m1);
    %ktt =u1+u2+u3+u4-d1-d2-d3-d4-d5
    %CG(1) = fact * exp((u1+u2+u3+u4-d1-d2-d3-d4-d5)/2);
    
    Sn(1)=fact; % sign;
    Kn(1) = (u1+u2+u3+u4-d1-d2-d3-d4-d5)/2; % exponents;
    CG(1)=Sn(1)*exp(Kn(1));
    
    %b=5
elseif (( abs(n1-n2)<abs(m3))&& (m3<0))
    % mp == -mp;
    % m == -m;
    u1 = logfact(-2*m3+1);
    u2 = logfact(n2+n1+m3);
    u3 = logfact(n1-m1);
    u4 = logfact(n2-m3+m1);
    d1 = logfact(n1+n2-m3+1);
    d2 = logfact(n1-n2-m3);
    d3 = logfact(n2-n1-m3);
    d4 = logfact(n1+m1);
    d5 = logfact(n2+m3-m1);
    fact = (-1)^(n1-m1);
    fact2 = (-1)^(n1+n2+m3);
    %ktt =u1+u2+u3+u4-d1-d2-d3-d4-d5
    %he=exp(ktt);
    %CG(1) = fact * fact2*exp((u1+u2+u3+u4-d1-d2-d3-d4-d5)/2);
    
    Sn(1)=fact * fact2; % sign;
    Kn(1) = (u1+u2+u3+u4-d1-d2-d3-d4-d5)/2; % exponents;
    CG(1)=Sn(1)*exp(Kn(1));
    %b=6
end

if (N3max<=10) % upward only; 
       % b=8
    k=1;
    ntop=N3max;
    for np= (N3min+1):ntop
        uA = 4*(np^2)*(2*np+1)*(2*np-1);
        dA = (np+m3)*(np-m3)*(n2-n1+np)*(n1-n2+np)*(n1+n2-np+1)*(n1+n2+np+1);
        A = sqrt(uA/dA);
        
        uB = (2*m1-m3)*np*(np-1) - (m3*n1)*(n1+1) + (m3*n2)*(n2+1);
        dB = 2*np*(np-1);
        B = uB/dB;
        
        uC = (np-m3-1)*(np+m3-1)*(n2-n1+np-1)*(n1-n2+np-1)*(n1+n2-np+2)...
            *(n1+n2+np);
        dC = 4.0*(np-1)^2*(2*np-3)*(2*np-1);
        C = sqrt(uC/dC);
        
        if (k==1)
            % value
            tp=A*B*Sn(1);
            % sign
            Sn(2)=sign(tp);
            % exponents
            Kn(2)=Kn(1)+log(abs(tp));
            CG(2)=Sn(2)*exp(Kn(2));
            k=k+1;
            %CG(2)=A*(B*CG(1));k= k+1;
            
        else
            
            if(Sn(k)~=0)&&(Sn(k-1)~=0)
                % value
                Dk=Kn(k)-Kn(k-1); % exponent difference;
                tp=A*(B*Sn(k)-C*Sn(k-1)*exp(-Dk));           
                Sn(k+1)=sign(tp);
                Kn(k+1)=Kn(k)+log(abs(tp));
                %CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                CG(k+1)=Sn(k+1)*exp(Kn(k+1));
                k = k+1;
                
            else
                
                CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                Kn(k+1)=log(abs(CG(k+1)));
                Sn(k+1)=sign(CG(k+1));
                k=k+1;
                                
            end
        end
        
    end
    %return;
end
if (N3max > 10) % combine upawrd and downward; 
        %b=9
        % let us do upward first;
        ntop=N3min+Nmid+1;
        k=1;
        %t=CG(1)
        for np = (N3min+1):ntop % upward
            
            uA = 4*(np^2)*(2*np+1)*(2*np-1);
            dA = (np+m3)*(np-m3)*(n2-n1+np)*(n1-n2+np)*(n1+n2-np+1)*(n1+n2+np+1);
            A = sqrt(uA/dA);
            
            uB = (2*m1-m3)*np*(np-1) - (m3*n1)*(n1+1) + (m3*n2)*(n2+1);
            dB = 2*np*(np-1);
            B = uB/dB;
            
            uC = (np-m3-1)*(np+m3-1)*(n2-n1+np-1)*(n1-n2+np-1)*(n1+n2-np+2)...
                *(n1+n2+np);
            dC = 4.0*(np-1)^2*(2*np-3)*(2*np-1);
            C = sqrt(uC/dC);
            
            if (k == 1)       
                % value
                tp=A*B*Sn(1);
                % sign
                Sn(2)=sign(tp);
                % exponents
                Kn(2)=Kn(1)+log(abs(tp));
                CG(2)=Sn(2)*exp(Kn(2));
                k=k+1;
                %  CG(2)=A*(B*CG(1));k= k+1;
            else
                               
                if(Sn(k-1)~=0)&&(Sn(k)~=0)  % advoid zeros;
                    % value
                    %k;
                    Dk=Kn(k)-Kn(k-1); % exponent difference;
                    tp=A*(B*Sn(k)-C*Sn(k-1)*exp(-Dk));
                    Sn(k+1)=sign(tp);
                    Kn(k+1)=Kn(k)+log(abs(tp));                    
                    %CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                    CG(k+1)=Sn(k+1)*exp(Kn(k+1));
                    k = k+1;
                                        
                else
                    %kkk=1
                    CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                    Sn(k+1) = sign(CG(k+1));
                    Kn(k+1) = log(abs(CG(k+1)));
                    k=k+1;    
                end
                % CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                % k = k+1;
            end
            
        end
        
        
        % start from N3max;
        k = nC;
        u1 = logfact(2*n1);
        u2 = logfact(2*n2);
        u3 = logfact(n2+n1+m3);
        u4 = logfact(n2+n1-m3);
        d1 = logfact(2*n1 + 2*n2);
        d2 = logfact(n1+m1);
        d3 = logfact(n1-m1);
        d4 = logfact(n2+m3-m1);
        d5 = logfact(n2-m3+m1);
        
        % sign
        Sn(k)=1.0;
        % exponent
        Kn(k)= (u1+u2+u3+u4-d1-d2-d3-d4-d5)/2;
        % CG
        CG(k) = Sn(k)*exp(Kn(k));
        
        %CG(k) = exp((u1+u2+u3+u4-d1-d2-d3-d4-d5)/2);
        
        ntop=ntop+1;
       
        for np = (N3max-1):-1:(ntop)
            
            uD = 4.0*(np+1)^2*(2*np+3)*(2*np+1);
            dD = (np+m3+1)*(np-m3+1)*(n2-n1+np+1)...
                *(n1-n2+np+1)*(n1+n2-np)*(n1+n2+np+2);
            D =sqrt(uD/dD);
            
            uE = (2*m1-m3)*(np+2)*(np+1)- (m3*n1)*(n1+1) +(m3*n2)*(n2+1);
            dE = 2.0*(np+2)*(np+1);
            E = uE/dE;
            
            uF = (np-m3+2)*(np+m3+2)*(n2-n1+np+2)*(n1-n2+np+2)*(n1+n2-np-1)*(n1+n2+np+3);
            dF = 4*(np+2)^2*(2*np+5)*(2*np+3);
            F =sqrt(uF/dF);
            
            if(k==nC)
                
                % value;
                tp=D*E*Sn(k);
                % sign
                Sn(k-1)=sign(tp);
                % exponents;
                Kn(k-1)=Kn(k)+log(abs(tp));
                % CGs
                CG(k-1)=Sn(k-1)*exp(Kn(k-1));
                k=k-1;
                
                %CG(k-1) = D*E*CG(k); k=k-1;
            else
                if(Sn(k+1)~=0)&&(Sn(k)~=0)  % advoid zeros;    
                    % values;
                    Dk=Kn(k)-Kn(k+1);
                    tp = D*(E*Sn(k)-F*Sn(k+1)*exp(-Dk));
                    % sign and exponents;
                    Sn(k-1)=sign(tp);
                    Kn(k-1)=Kn(k)+log(abs(tp));
                    CG(k-1)=Sn(k-1)*exp(Kn(k-1));
                    k=k-1;
                else
                    CG(k-1) = D*(E*CG(k) -F*CG(k+1));
                    Kn(k-1)=log(abs(CG(k-1)));
                    Sn(k-1)=sign(CG(k-1));
                    k=k-1;
                                        
                end
                   %CG(k-1) = D*(E*CG(k) -F*CG(k+1)); k= k-1;
            end
            
        end
        
end


for ik = n3a:n3b
    if (ik<N3min)
        cg(ik-n3a+1)=0.0;
    elseif (ik>N3max)
        cg(ik-n3a+1)=0.0;
    else
        %cgn(ik-m2_a+1) = Cx(ik-m2_min+1);
        %m3=m1+ik;
        cg(ik-n3a+1)=CG(ik-N3min+1); 
        snt(ik-n3a+1)=Sn(ik-N3min+1);
        knt(ik-n3a+1)=Kn(ik-N3min+1);
        %cgn(ik-m2_a+1) = Cx(ik-m2_min+1)*sqrt(2*j3+1)*(-1)^(j1-j2+m3);
    end
    
end

%b=10;
end