t = 0.00001; %s
%T = 0.1; %s
Cd = 0.6;
D = 0.01; %m
H2 = 1; %g
AF = 1;
TR = 293; %k
Tref = 298; %k
lengthcc = 0.25; %m
diametercc = 0.07; %m
VT = 0.06; %m3
Tt0 = 293; %k
Pt0 = 101325; %pa
realgas = 1;

aH2O = 0.5537; %Pa.m6/mol2
aO2 = 0.1382; %Pa.m6/mol2
aN2 = 0.137; %Pa.m6/mol2
bH2O = 0.00003049; %m3/mol
bO2 = 0.00003186; %m3/mol
bN2 = 0.0000387; %m3/mol

Ru = 8.31447; %kj/kmol.K

MMH2 = 2.01588; %kg/kmol
MMO2 = 31.999; %kg/kmol
MMN2 = 28.013; %kg/kmol
MMH2O = 18.015; %kg/kmol
MMair = 28.97; %kg/kmol

hfH2O = -241820; %kj/mol

H2OCpT298 = 1.8723; %kj/kg.k
O2CpT298 = 0.918; %kj/kg.k
N2CpT298 = 1.039; %kj/kg.k

H2CpTR = 14.4; %kj/kg.k
O2CpTR = 0.917; %kj/kg.k
N2CpTR = 1.039; %kj/kg.k

nH2R = H2/MMH2; %mol
nO2R = (nH2R/2)*AF; %mol
nN2R = (nH2R/2)*AF*3.76; %mol
nH2OP = nH2R; %mol
nO2P = (nH2R/2)*(AF-1); %mol
nN2P = nN2R; %mol
VC = lengthcc*pi*((diametercc/2)^2); %m3

UR = nH2R*(H2CpTR*MMH2*(TR-Tref)-Ru*TR)+nO2R*(O2CpTR*MMO2*(TR-Tref)-Ru*TR)+nN2R*(N2CpTR*MMN2*(TR-Tref)-Ru*TR); %J

Tcc0 = TR;
CpH2O = H2OCpT298;
CpO2 = O2CpT298;
CpN2 = N2CpT298;
for j = 1:4
    Tcc0 = (UR-nH2OP*hfH2O - nH2OP*CpH2O*MMH2O*(-Tref) - nO2P*CpO2*MMO2*(-Tref) - nN2P*CpN2*MMN2*(-Tref))/(nH2OP*CpH2O*MMH2O - nH2OP*Ru + nO2P*CpO2*MMO2 - nO2P*Ru + nN2P*CpN2*MMN2 - nN2P*Ru); %K
    Tavg = (Tcc0+Tref)/2;
    CpH2O = (32.24+0.1923*(10^-2)*Tavg+1.055*(10^-5)*(Tavg^2)-3.595*(10^-9)*(Tavg^3))/MMH2O; %Kj/Kg.K
    CpO2 = (25.48+1.52*(10^-2)*Tavg-0.7155*(10^-5)*(Tavg^2)+1.312*(10^-9)*(Tavg^3))/MMO2; %Kj/Kg.K
    CpN2 = (28.9-0.1571*(10^-2)*Tavg+0.8081*(10^-5)*(Tavg^2)-2.873*(10^-9)*(Tavg^3))/MMN2; %Kj/Kg.K
end

if realgas == 0
    Pcc0 = ((nH2OP+nO2P+nN2P)*Ru*Tcc0)/VC; %Pa
else
    Pcc0 = (((nH2OP*Ru*Tcc0)/(VC-bH2O*nH2OP))-(aH2O*(nH2OP^2)/(VC^2)))+(((nO2P*Ru*Tcc0)/(VC-bO2*nO2P))-(aO2*(nO2P^2))/(VC^2))+((nN2P*Ru*Tcc0)/(VC-bN2*nN2P))-(aN2*nN2P^2)/(VC^2); %Pa
end

mH2O = nH2OP * MMH2O; %g
mO2 = nO2P * MMO2; %g
mN2 = nN2P * MMN2; %g

mtotal = mH2O + mO2 + mN2; %g
ntotal = nH2OP + nO2P + nN2P; %mol

MMM = mtotal / ntotal; %kg/kmol
R = Ru/MMM; %kj/kg.k
A = ((D/2)^2)*pi; %m2
rho = (mtotal*0.001)/VC; %kg/m3

YH2O = nH2OP/ntotal;
YO2 = nO2P/ntotal;
YN2 = nN2P/ntotal;
XH2O = mH2O/mtotal;
XO2 = mO2/mtotal;
XN2 = mN2/mtotal;

CpH2O = ((32.24+0.1923*(10^-2)*Tcc0+1.055*(10^-5)*(Tcc0^2)-3.595*(10^-9)*(Tcc0^3))/MMH2O);
CpO2 = ((25.48+1.52*(10^-2)*Tcc0-0.7155*(10^-5)*(Tcc0^2)+1.312*(10^-9)*(Tcc0^3))/MMO2);
CpN2 = ((28.9-0.1571*(10^-2)*Tcc0+0.8081*(10^-5)*(Tcc0^2)-2.873*(10^-9)*(Tcc0^3))/MMN2);
CvH2O = CpH2O-(Ru/MMH2O);
CvO2 = CpO2-(Ru/MMO2);
CvN2 = CpN2-(Ru/MMN2);
kH2O = CpH2O/CvH2O;
kO2 = CpO2/CvO2;
kN2 = CpN2/CvN2;
k = (((kH2O-1)*(kO2-1)*(kN2-1)*(ntotal))/((nH2OP*(kO2-1)*(kN2-1))+(nO2P*(kH2O-1)*(kN2-1))+(nN2P*(kH2O-1)*(kO2-1))))+1;

m0 = A*Cd*((k*Pcc0*rho)^0.5)*((2/(k+1))^((k+1)/(2*(k-1))))*1000; %g

M = zeros ( 10000,31);
M(1,1)= 0;
M(1,2)= m0;
M(1,3)= 0;
M(1,4)= Pt0*VT/((Ru/MMair)*Tt0);
M(1,5)= Tt0;
M(1,6)= Pt0;
M(1,7)= mtotal;
M(1,8)= Tcc0;
M(1,9)= Pcc0;
M(1,10)= ((32.24+0.1923*(10^-2)*M(1,5)+1.055*(10^-5)*(M(1,5)^2)-3.595*(10^-9)*(M(1,5)^3))/MMH2O);
M(1,11)= ((25.48+1.52*(10^-2)*M(1,5)-0.7155*(10^-5)*(M(1,5)^2)+1.312*(10^-9)*(M(1,5)^3))/MMO2);
M(1,12)= ((28.9-0.1571*(10^-2)*M(1,5)+0.8081*(10^-5)*(M(1,5)^2)-2.873*(10^-9)*(M(1,5)^3))/MMN2);
M(1,13)= M(1,10)-(Ru/MMH2O);
M(1,14)= M(1,11)-(Ru/MMO2);
M(1,15)= M(1,12)-(Ru/MMN2);
M(1,16)= ((32.24+0.1923*(10^-2)*M(1,8)+1.055*(10^-5)*(M(1,8)^2)-3.595*(10^-9)*(M(1,8)^3))/MMH2O);
M(1,17)= ((25.48+1.52*(10^-2)*M(1,8)-0.7155*(10^-5)*(M(1,8)^2)+1.312*(10^-9)*(M(1,8)^3))/MMO2);
M(1,18)= ((28.9-0.1571*(10^-2)*M(1,8)+0.8081*(10^-5)*(M(1,8)^2)-2.873*(10^-9)*(M(1,8)^3))/MMN2);
M(1,19)= M(1,16)-(Ru/MMH2O);
M(1,20)= M(1,17)-(Ru/MMO2);
M(1,21)= M(1,18)-(Ru/MMN2);
M(1,22)= M(1,3);

Y2H2O = (YH2O*M(1,22)/MMM)/(M(1,4)/MMair+M(1,22)/MMM);
Y2O2 = (((1/4.76)*M(1,4)/MMair)+(YO2*M(1,22)/MMM))/((M(1,4)/MMair)+(M(1,22)/MMM));
Y2N2 = (((3.76/4.76)*M(1,4)/MMair)+(YN2*M(1,22))/MMM)/((M(1,4)/MMair)+(M(1,22)/MMM));
X2H2O = (XH2O*M(1,22))/(M(1,4)+M(1,22));
X2O2 = ((((1/4.76)*M(1,4)/MMair)*MMO2)+ XO2*M(1,22))/(M(1,4)+M(1,22));
X2N2 = ((((3.76/4.76)*M(1,4)/MMair)*MMN2)+XN2*M(1,22))/(M(1,4)+M(1,22));

M(1,23)= Y2H2O;
M(1,24)= Y2O2;
M(1,25)= Y2N2;
M(1,26)= X2H2O;
M(1,27)= X2O2;
M(1,28)= X2N2;

for i = 2:100000
    M(i,1)= M(i-1,1)+t ; 
    
    M(i,3) = M(i-1,2)*t;
    M(i,22) = M(i-1,22)+M(i,3);
    
    M(i,4) = M(i-1,4)+ M(i,3);
    
    Y2H2O = (YH2O*M(i,22)/MMM)/(M(1,4)/MMair+M(i,22)/MMM);
    Y2O2 = (((1/4.76)*M(1,4)/MMair)+(YO2*M(i,22)/MMM))/((M(1,4)/MMair)+(M(i,22)/MMM));
    Y2N2 = (((3.76/4.76)*M(1,4)/MMair)+(YN2*M(i,22))/MMM)/((M(1,4)/MMair)+(M(i,22)/MMM));
    X2H2O = (XH2O*M(i,22))/(M(1,4)+M(i,22));
    X2O2 = ((((1/4.76)*M(1,4)/MMair)*MMO2)+ XO2*M(i,22))/(M(1,4)+M(i,22));
    X2N2 = ((((3.76/4.76)*M(1,4)/MMair)*MMN2)+XN2*M(i,22))/(M(1,4)+M(i,22));
    
    M(i,23)= Y2H2O;
    M(i,24)= Y2O2;
    M(i,25)= Y2N2;
    M(i,26)= X2H2O;
    M(i,27)= X2O2;
    M(i,28)= X2N2;
    
    Tavg1 = (M(i-1,8)+Tref)/2;
    CpH2Oavg1 = (32.24+0.1923*(10^-2)*Tavg1+1.055*(10^-5)*(Tavg1^2)-3.595*(10^-9)*(Tavg1^3))/MMH2O; %Kj/Kg.K
    CpO2avg1 = (25.48+1.52*(10^-2)*Tavg1-0.7155*(10^-5)*(Tavg1^2)+1.312*(10^-9)*(Tavg1^3))/MMO2; %Kj/Kg.K
    CpN2avg1 = (28.9-0.1571*(10^-2)*Tavg1+0.8081*(10^-5)*(Tavg1^2)-2.873*(10^-9)*(Tavg1^3))/MMN2; %Kj/Kg.K

    Tavg2 = (M(i-1,5)+Tref)/2;
    CpH2Oavg2 = (32.24+0.1923*(10^-2)*Tavg1+1.055*(10^-5)*(Tavg1^2)-3.595*(10^-9)*(Tavg1^3))/MMH2O; %Kj/Kg.K
    CpO2avg2 = (25.48+1.52*(10^-2)*Tavg1-0.7155*(10^-5)*(Tavg1^2)+1.312*(10^-9)*(Tavg1^3))/MMO2; %Kj/Kg.K
    CpN2avg2 = (28.9-0.1571*(10^-2)*Tavg1+0.8081*(10^-5)*(Tavg1^2)-2.873*(10^-9)*(Tavg1^3))/MMN2; %Kj/Kg.K

    M(i,5) = (M(i,3)*(XH2O*(hfH2O+CpH2Oavg1*(M(i-1,8)-Tref))+XO2*CpO2avg1*(M(i-1,8)-Tref)+XN2*CpN2avg1*(M(i-1,8)-Tref))+M(i-1,4)*(M(i-1,28)*(CpN2avg2*(M(i-1,5)-Tref)-(Ru/MMN2)*M(i-1,5))+M(i-1,27)*(CpO2avg2*(M(i-1,5)-Tref)-(Ru/MMO2)*M(i-1,5))+M(i-1,26)*(hfH2O+CpH2Oavg2*(M(i-1,5)-Tref)-(Ru/MMH2O)*M(i-1,5)))+M(i,28)*M(i,4)*CpN2avg2*Tref+M(i,27)*M(i,4)*CpO2avg2*Tref+M(i,26)*M(i,4)*CpH2Oavg2*Tref-M(i,26)*M(i,4)*hfH2O)/(M(i,4)*(M(i,28)*(CpN2avg2-(Ru/MMN2))+M(i,27)*(CpO2avg2-(Ru/MMO2))+M(i,26)*(CpH2Oavg2-(Ru/MMH2O))));

    if realgas == 0
        M(i,6) = M(i,4)*R*M(i,5)/VT;
    else
        M(i,6) = (((((M(i,22)*YH2O)/MMM)*Ru*M(i,5))/(VT-bH2O*((M(i,22)*YH2O)/MMM)))-(aH2O*(((M(i,22)*YH2O)/MMM)^2)/(VT^2)))+((((((M(i,22)*YO2)/MMM)+((1/4.76)*M(1,4)/MMair))*Ru*M(i,5))/(VT-bO2*(((M(i,22)*YO2)/MMM)+((1/4.76)*M(1,4)/MMair))))-(aO2*((((M(i,22)*YO2)/MMM)+((1/4.76)*M(1,4)/MMair))^2)/(VT^2)))+((((((M(i,22)*YN2)/MMM)+((3.76/4.76)*M(1,4)/MMair))*Ru*M(i,5))/(VT-bN2*(((M(i,22)*YN2)/MMM)+((3.76/4.76)*M(1,4)/MMair))))-(aN2*((((M(i,22)*YN2)/MMM)+((3.76/4.76)*M(1,4)/MMair))^2)/(VT^2))) ;
    end

    M(i,7)= M(i-1,7)-M(i,3);
    
    M(i,8)= (-M(i,3)*(XH2O*(hfH2O+CpH2Oavg1*(M(i-1,8)-Tref))+XO2*CpO2avg1*(M(i-1,8)-Tref)+XN2*CpN2avg1*(M(i-1,8)-Tref)) + M(i-1,7)*(XN2*(CpN2avg1*(M(i-1,8)-Tref)-(Ru/MMN2)*M(i-1,8))+XO2*(CpO2avg1*(M(i-1,8)-Tref)-(Ru/MMO2)*M(i-1,8))+XH2O*(hfH2O + CpH2Oavg1*(M(i-1,8)-Tref)-(Ru/MMH2O)*M(i-1,8)))+XN2*M(i,7)*CpN2avg1*Tref+XO2*M(i,7)*CpO2avg1*Tref+XH2O*M(i,7)*CpH2Oavg1*Tref-XH2O*M(i,7)*hfH2O)/(M(i,7)*(XN2*(CpN2avg1-(Ru/MMN2))+XO2*(CpO2avg1-(Ru/MMO2))+XH2O*(CpH2Oavg1-(Ru/MMH2O))));
    
    if realgas == 0
        M(i,9)= M(i,7)*R*1000*M(i,8)/VC;
    else
        M(i,9)= (((((M(i,7)*YH2O)/MMM)*Ru*M(i,8))/(VC-bH2O*((M(i,7)*YH2O)/MMM)))-(aH2O*(((M(i,7)*YH2O)/MMM)^2)/(VC^2)))+(((((M(i,7)*YO2)/MMM)*Ru*M(i,8))/(VC-bO2*((M(i,7)*YO2)/MMM)))-(aO2*(((M(i,7)*YO2)/MMM)^2)/(VC^2)))+(((((M(i,7)*YN2)/MMM)*Ru*M(i,8))/(VC-bN2*((M(i,7)*YN2)/MMM)))-(aN2*(((M(i,7)*YN2)/MMM)^2)/(VC^2)));
    end
    
    M(i,10)= ((32.24+0.1923*(10^-2)*M(i,5)+1.055*(10^-5)*(M(i,5)^2)-3.595*(10^-9)*(M(i,5)^3))/MMH2O);
    M(i,11)= ((25.48+1.52*(10^-2)*M(i,5)-0.7155*(10^-5)*(M(i,5)^2)+1.312*(10^-9)*(M(i,5)^3))/MMO2);
    M(i,12)= ((28.9-0.1571*(10^-2)*M(i,5)+0.8081*(10^-5)*(M(i,5)^2)-2.873*(10^-9)*(M(i,5)^3))/MMN2);
    M(i,13)= M(i,10)-(Ru/MMH2O);
    M(i,14)= M(i,11)-(Ru/MMO2);
    M(i,15)= M(i,12)-(Ru/MMN2);
    M(i,16)= ((32.24+0.1923*(10^-2)*M(i,8)+1.055*(10^-5)*(M(i,8)^2)-3.595*(10^-9)*(M(i,8)^3))/MMH2O);
    M(i,17)= ((25.48+1.52*(10^-2)*M(i,8)-0.7155*(10^-5)*(M(i,8)^2)+1.312*(10^-9)*(M(i,8)^3))/MMO2);
    M(i,18)= ((28.9-0.1571*(10^-2)*M(i,8)+0.8081*(10^-5)*(M(i,8)^2)-2.873*(10^-9)*(M(i,8)^3))/MMN2);
    M(i,19)= M(i,16)-(Ru/MMH2O);
    M(i,20)= M(i,17)-(Ru/MMO2);
    M(i,21)= M(i,18)-(Ru/MMN2);
    
    kH2O = M(i,16)/M(i,19);
    kO2 = M(i,17)/M(i,20);
    kN2 = M(i,18)/M(i,21);
    %k = (((kH2O-1)*(kO2-1)*(kN2-1)*(M(i,7/MMM))/(((YH2O*M(i,7)/MMM)*(kO2-1)*(kN2-1))+((YO2*M(i,7)/MMM)*(kH2O-1)*(kN2-1))+((YN2*M(i,7)/MMM)*(kH2O-1)*(kO2-1)))))+1;
    %k = (((kH2O-1)*(kO2-1)*(kN2-1)*(nH2OP+nO2P+nN2P))/((nH2OP*(kO2-1)*(kN2-1))+(nO2P*(kH2O-1)*(kN2-1))+(nN2P*(kH2O-1)*(kO2-1))))+1;
    k = (((kH2O-1)*(kO2-1)*(kN2-1)*(M(i,7)/MMM))/(YH2O*(M(i,7)/MMM)*(kO2-1)*(kN2-1)+YO2*(M(i,7)/MMM)*(kH2O-1)*(kN2-1)+YN2*(M(i,7)/MMM)*(kO2-1)*(kH2O-1)))+1;
    M(i,32) = k;
    
    M(i,29) = M(i,6)/M(i,9);
    M(i,30) = ((2/(k+1))^(k/(k-1)));
    rho = (M(i,7)/1000)/VC;
    
    %M(i,2)= A*Cd*((k*M(i,9)*rho)^0.5)*((2/(k+1))^((k+1)/(2*(k-1))))*1000;
    if M(i,29) <= M(i,30)
    M(i,2)= A*Cd*((k*M(i,9)*rho)^0.5)*((2/(k+1))^((k+1)/(2*(k-1))))*1000;
    %M(i,2)= A*Cd*((k*M(i,9)*rho*((2/(k+1))^((k+1)/(k-1))))^0.5)*1000;
    M(i,31) = 1;
    else 
    M(i,2)= A*Cd*((2*M(i,9)*rho*(k/(k-1))*(((M(i,6)/M(i,9))^(2/k))-((M(i,6)/M(i,9))^((k+1)/k))))^0.5);
    M(i,31) = 0;
    end




end
