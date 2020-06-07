function  [cg, snt,  knt] = CG_coefficients_n3_log(n1, n2, n3, m1, m2)

% Code Purpose:
% Compute particular Clebsch-Gordan coefficient using
% the sign-exponent recurrence method;;

% Inputs:
% n1, n2, n3; -------------% principle quantum numbers 
% m1  -----------------% magnetic quantum numbers for j1
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

if (n3>N3max)||(n3<N3min)
    cg=0.0;
    fprintf('j3 is not within the range.');
    return;
end

if (n1<0)||(n2<0)
    
    cg=0.0;
    fprintf('n1 and n2 shoud be positive number.');
    return;
end
if (abs(m2)>n2) || (abs(m1)>n1)
    fprintf('m1 and m2 are not within the range.');
    cg = 0.0; return;
end

if(N3max==0)
    fprintf('n3 max is zero.');
    cg=0.0; return;
end

nC = N3max-N3min+1;
Nmid = round((N3min + N3max)/2.);
CG(nC) = zeros; % CG-coefficients;
Sn(nC) = ones;  % sign of coefficients;
Kn(nC) = zeros; % exponents of the coefficients;

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

if (nC==1)||(n3==N3min)
    cg = Sn(1)*exp(Kn(1));
    snt=Sn(1);
    knt=Kn(1);
    %b=7
    return;
end
%t=CG(1);

if (N3max<=30)
       % b=8
    k=1;
    ntop=n3;
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
    cg=CG(ntop-N3min+1);
    snt=Sn(ntop-N3min+1);
    knt=Kn(ntop-N3min+1);
    return;
end
if (N3max > 30)
        %b=9
    if (n3<=Nmid) % upward;
        ntop=n3;
        k=1;
        %t=CG(1)
        for np = (N3min+1):ntop
            
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
                    k;
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
                    Kn(k+1) = exp(abs(CG(k+1)));
                    k=k+1;    
                end
                % CG(k+1) = A*(B*CG(k) - C*CG(k-1));
                % k = k+1;
            end
            
        end
        
        cg=CG(ntop-N3min+1);
        snt=Sn(ntop-N3min+1);
        knt=Kn(ntop-N3min+1);
        %nn=10;
        %plot(1:length(Kn), Kn); 
        return;
        
    else        % downward;
        
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
        
        ntop=n3;
        
        if(ntop==N3max)
            cg= CG(k);
            snt=Sn(k);
            knt=Kn(k);
            return;
        end
        
        
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
        
        cg=CG(ntop-N3min+1);
        snt=Sn(ntop-N3min+1);
        knt=Kn(ntop-N3min+1);
        return;
        
    end
    
end

%b=10;
end