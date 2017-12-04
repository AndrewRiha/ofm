% Filename: SRT4_Propulsion_OxidizerFill
% Created By: Andrew Riha
% Last Edited By: Andrew Riha
% Last Edited On: November 3, 2017
% Project: Oxidizer Fill
% Description: This code is the full model for filling the run tank with
% oxidizer. Oxidizer code models the run tank as initially full of Nitrous
% oxide at ambient temperature and pressure. Assumptions include: all
% boiling happens instantaneously, the exit plane of the vent line has the
% smallest area in the vent line, all parts of each tank instantaneously
% feel their respective temperature

% The equations and constants used come
% from "Fundamentals of Engineering Thermodynamics" Chapter 11.4, "Rocket 
% Propulsion Elements" by Sutton Chapter 8.1, "Thermophysical properties
% of nitrous oxide (ESDU 91022)", NIST's webpage on Nitrous oxide data 
% "http://webbook.nist.gov/cgi/cbook.cgi?ID=C10024972&Mask=1", The
% Engineering ToolBox's webpages on gases
% "http://www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html"
% and dry air
% "http://www.engineeringtoolbox.com/molecular-mass-air-d_679.html"



% Initial Conditions
    % "1" variables refer to the fill tank, "2" variables refer to run 
    % tank, "12" variables refer to the transfer between the fill tank and
    % run tank, "23" variables refer to the transfer between the run tank
    % and the atmosphere
clear; clc; close all
TAmb=70; % Ambient Temperature [deg F]
T1=75; % Fill tank temperature
pAmb=1; % Ambient Pressure [atm]
V1=1.553845; % Volume [ft^3]
w1tot=50; % Weight of contents of fill tank [lbf]
V2=0.318287; % Volume [ft^3]
w2Target=10; % Desired weight of nitrous oxide in run tank [lbf]
Cd12=1; % Discharge coefficient of hose connecting tanks PLACEHOLDER
A12=(1/4)^2*pi/4; % Inside area of hose connecting tanks [in^2] PLACEHOLDER
A23=(1/16)^2*pi/4; % Inside area of hose connecting run tank to atmosphere
Tc=557.23; % Critical temperature of Nitrous oxide [deg R]
pc=151400; % Critical pressure of Nitrous oxide [psf]
rhoc=0.877; % Critical density of Nitrous oxide [slug/ft^3]
RUniv=3.40690602232; % Universal Gas Constant [ft*lbf/((deg R)*mol)]
molecMassN2O=3.01584847e-3; % Molecular Mass of Nitrous oxide [slug/mol]
gammaN2O=1.27; % Specific Heat Ratio of Nitrous oxide
gammaN2OKT=4/3; % Specific Heat Ratio of Nitrous oxide from Kinetic Theory


% Read in Thermodynamic Data for N2O
load N2OThermProp
% p: Vapor pressure [lbf/ft^2]
% T: Temperature [deg R]
% nuf: Saturated liquid specific volume [ft^3/mol]
% nug: Saturated vapor specific volume [ft^3/mol]
% hf: Saturated liquid enthalpy [ft*lbf/mol]
% hg: Saturated vapor enthalpy [ft*lbf/mol]
% sf: Saturated liquid entropy [ft*lbf/(mol*deg R)]
% sg: Saturated vapor entropy [ft*lbf/(mol*deg R)]
% sfRef: integral(cpf/T*dT,TRef,T(i)), where TRef=328.194 deg R, and T(i) 
    % is the temperature at which sfRef is being calculated
    % [ft*lbf/(mol*deg R)]
% sgRef: integral(cpf/T*dT,TRef,T(i)), where TRef=328.194 deg R, and T(i) 
    % is the temperature at which sgRef is being calculated
    % [ft*lbf/(mol*deg R)]
% cpf: Saturated liquid constant pressure specific heat capacity
    % [ft*lbf/(mol*deg R)]
% cpg: Saturated vapor constant pressure specific heat capacity
    % [ft*lbf/(mol*deg R)]
% TGas: Temperature [deg R]
% pGas: Pressure [lbf/ft^2]
% sGas: Table of entropy for a given temperature and pressure. Temperature
    % is given by TGas and corresponds with the rows of sGas. Pressure
    % is given by pGas and corresponds with the columns of sGas.
    % [ft*lbf/(mol*deg R)]
%***ADD sGasRef FROM NIST TABLES***

    
% Initial Conditions
i=1; % Counter
tStep=1e-2; % [s]
TStep=.01; % [deg R]

TAmb=TAmb+459.67; % [deg R]
T1=T1+459.67; % Temperature [deg R]
T2=TAmb; % Temperature [deg R]

% Run Tank
iInterp2=find(T>T2,1); % First index past current temperature
interp2=(T(iInterp2)-T2)/(T(iInterp2)-T(iInterp2-1)); % Linear 
    % interpolation
    
m2Target=w2Target/32.174; % Desired mass of Nitrous oxide in run tank 
    % [slug]
w2tot=2.6; % [lbf] ***PLACEHOLDER***
m2tot=w2tot/32.174; % Total mass of Nitrous oxide in run tank [slug]
n2tot=m2tot/molecMassN2O; % [mol]

p2=(p(iInterp2)-p(iInterp2-1))*interp2+p(iInterp2-1); % [psf]

nu2f=(nuf(iInterp2)-nuf(iInterp2-1))*interp2+nuf(iInterp2-1); % [ft^3/mol]
nu2g=(nug(iInterp2)-nug(iInterp2-1))*interp2+nug(iInterp2-1); % [ft^3/mol]
nu2=V2/m2tot*molecMassN2O; % Specific volume [ft^3/slug]
x2=(nu2-nu2f)/(nu2g-nu2f); % Quality in fill tank

pAmb=pAmb*2116.22; % [psf]

n2g=x2*n2tot; % [mol] Moles of gas in fill tank
n2f=n2tot-n2g; % [mol]  Moles of liquid in fill tank

h2g(i)=(hg(iInterp2)-hg(iInterp2-1))*interp2+hg(iInterp2-1); % ***PLACEHOLDER***

cp2f(i)=(cpf(iInterp2)-cpf(iInterp2-1))*interp2+cpf(iInterp2-1);
cp2g(i)=(cpg(iInterp2)-cpg(iInterp2-1))*interp2+cpg(iInterp2-1);

nu2f(i)=(nuf(iInterp2)-nuf(iInterp2-1))*interp2+cpf(iInterp2-1);
nu2g(i)=(nug(iInterp2)-nug(iInterp2-1))*interp2+nug(iInterp2-1);

iInterp2d1=find(T2(i)-TStep<T,1); 
interp2d1=(T(iInterp2d1)-(T2(i)-TStep))/...
    (T(iInterp2d1)-T(iInterp2d1-1));
iInterp2d2=find(T2(i)+TStep<T,1); 
interp2d2=(T(iInterp2d2)-(T2(i)+TStep))/...
    (T(iInterp2d2)-T(iInterp2d2-1));

p2d1=(p(iInterp2d1)-p(iInterp2d1-1))*interp2d1+p(iInterp2d1-1);
p2d2=(p(iInterp2d2)-p(iInterp2d2-1))*interp2d2+p(iInterp2d2-1);

cp2fd1=(cpf(iInterp2d1)-cpf(iInterp2d1-1))*interp2d1+cpf(iInterp2d1-1);
cp2fd2=(cpf(iInterp2d2)-cpf(iInterp2d2-1))*interp2d2+cpf(iInterp2d2-1);
cp2gd1=(cpg(iInterp2d1)-cpg(iInterp2d1-1))*interp2d1+cpg(iInterp2d1-1);
cp2gd2=(cpg(iInterp2d2)-cpg(iInterp2d2-1))*interp2d2+cpg(iInterp2d2-1);

nu2gd1=(nug(iInterp2d1)-nug(iInterp2d1-1))*interp2d1+nug(iInterp2d1-1);
nu2gd2=(nug(iInterp2d2)-nug(iInterp2d2-1))*interp2d2+nug(iInterp2d2-1);

dp2dT(i)=(p2d2-p2d1)/(2*TStep);

dcp2fdT(i)=(cp2fd2-cp2fd1)/(2*TStep);
dcp2gdT(i)=(cp2gd2-cp2gd1)/(2*TStep);

dnu2gdT(i)=(nu2gd2-nu2gd1)/(2*TStep);


% Fill Tank
iInterp1=find(T>T1,1); % First index past current temperature
interp1=(T(iInterp1)-T1)/(T(iInterp1)-T(iInterp1-1)); % Linear 
    % interpolation

m1tot=w1tot/32.174; % Total mass of Nitrous oxide in fill tank [slug]
n1tot=m1tot/molecMassN2O; % [mol]

p1(i)=(p(iInterp1)-p(iInterp1-1))*interp1+p(iInterp1-1); % [psf]

nu1f(i)=(nuf(iInterp1)-nuf(iInterp1-1))*interp1+nuf(iInterp1-1); % [ft^3/mol]
nu1g(i)=(nug(iInterp1)-nug(iInterp1-1))*interp1+nug(iInterp1-1); % [ft^3/mol]
nu1(i)=V1/m1tot*molecMassN2O; % Specific volume [ft^3/slug]
x1(i)=(nu1(i)-nu1f(i))/(nu1g(i)-nu1f(i)); % Quality in fill tank

n1g(i)=x1(i)*n1tot; % [mol] Moles of gas in fill tank
n1f(i)=n1tot-n1g(i); % [mol]  Moles of liquid in fill tank

h1f(i)=(hf(iInterp1)-hf(iInterp1-1))*interp1+hf(iInterp1-1); % ***PLACEHOLDER***

cp1f(i)=(cpf(iInterp1)-cpf(iInterp1-1))*interp1+cpf(iInterp1-1);
cp1g(i)=(cpg(iInterp1)-cpg(iInterp1-1))*interp1+cpg(iInterp1-1);

iInterp1d1=find(T1(i)-TStep<T,1); 
interp1d1=(T(iInterp1d1)-(T1(i)-TStep))/...
    (T(iInterp1d1)-T(iInterp1d1-1));
iInterp1d2=find(T1(i)+TStep<T,1); 
interp1d2=(T(iInterp1d2)-(T1(i)+TStep))/...
    (T(iInterp1d2)-T(iInterp1d2-1));

p1d1=(p(iInterp1d1)-p(iInterp1d1-1))*interp1d1+p(iInterp1d1-1);
p1d2=(p(iInterp1d2)-p(iInterp1d2-1))*interp1d2+p(iInterp1d2-1);

cp1fd1=(cpf(iInterp1d1)-cpf(iInterp1d1-1))*interp1d1+cpf(iInterp1d1-1);
cp1fd2=(cpf(iInterp1d2)-cpf(iInterp1d2-1))*interp1d2+cpf(iInterp1d2-1);
cp1gd1=(cpg(iInterp1d1)-cpg(iInterp1d1-1))*interp1d1+cpg(iInterp1d1-1);
cp1gd2=(cpg(iInterp1d2)-cpg(iInterp1d2-1))*interp1d2+cpg(iInterp1d2-1);

nu1gd1=(nug(iInterp1d1)-nug(iInterp1d1-1))*interp1d1+nug(iInterp1d1-1);
nu1gd2=(nug(iInterp1d2)-nug(iInterp1d2-1))*interp1d2+nug(iInterp1d2-1);

dp1dT(i)=(p1d2-p1d1)/(2*TStep);

dcp1fdT(i)=(cp1fd2-cp1fd1)/(2*TStep);
dcp1gdT(i)=(cp1gd2-cp1gd1)/(2*TStep);

dnu1gdT(i)=(nu1gd2-nu1gd1)/(2*TStep);

% Transient Initial Conditions
t=0; % [s]
A12=A12/144; % [ft^2]
A23=A23/144; % [ft^2]
mdot12=A12*Cd12/nu1f(i)*sqrt(2*(p1(i)-p2(i))*nu1f(i)); % [slug/s]
mdot23=mdot12*0.05; % [slug/s] ***PLACEHOLDER***
eps=10e-3; % Desired accuracy
err=1; % Error in temperature
checkSat=false; % Saturation Check
regime=0; % Choked flow marker. If regime==1, vent line flow is choked, if 
    % regime==0, vent line flow is subsonic

% Transfer mass    
mTransfer12(i)=mdot12(i)*tStep; % [slug]
mTransfer23(i)=mdot23(i)*tStep; % [slug]
m1tot(i)=m1tot(i)-mTransfer12(i); % [slug]
m2tot=m2tot+mTransfer12(i)-mTransfer23(i); % [slug]

ndot12(i)=mdot12(i)/molecMassN2O; % [mol/s]
ndot23(i)=mdot23(i)/molecMassN2O; % [mol/s]

nTransfer12(i)=mTransfer12(i)/molecMassN2O; % [mol]
nTransfer23(i)=mTransfer23(i)/molecMassN2O; % [mol]
    
% Introduce Nitrous oxide to Run Tank
while m2tot<m2Target
    % Compute Derivatives in time
        % A will have three rows related to the mass equation, the energy 
            % equation, and the derivative of the ideal gas law. The three
            % columns will represent the coefficients on dnfdt, dngdt, and
            % dTdt. B is a vector that represents the constants with
            % repsect to these variables.
    A2(1,:)=[1 1 0];
    
    A2(2,1)=cp2f(i)*T2(i);
    A2(2,2)=cp2g(i)*T2(i);
    A2(2,3)=n2f(i)*(dcp2fdT(i)*T2(i)+cp2f(i))+...
        n2g(i)*(dcp2gdT(i)*T2(i)+cp2g(i));
    
    A2(3,1)=0;
    A2(3,2)=p2(i)*nu2g(i)+RUniv*T2(i);
    A2(3,3)=dp2dT(i)*(n2g(i)*nu2g(i)-V2)+p2(i)*n2g(i)*dnu2gdT(i)+...
        n2g(i)*RUniv;
    
    B2=[ndot12(i)-ndot23(i); (ndot12(i)*h1f(i)-ndot23(i)*h2g(i))*gammaN2O; 0];
    
    derivT2=A2\B2;
    
    A2Check(i)=cond(A2);
    
    
    A1(1,:)=[1 1 0];
    
    A1(2,1)=cp1f(i)*T1(i);
    A1(2,2)=cp1g(i)*T1(i);
    A1(2,3)=n1f(i)*(dcp1fdT(i)*T1(i)+cp1f(i))+...
        n1g(i)*(dcp1gdT(i)*T1(i)+cp1g(i));
    
    A1(3,1)=0;
    A1(3,2)=p1(i)*nu1g(i)+RUniv*T1(i);
    A1(3,3)=dp1dT(i)*(n1g(i)*nu1g(i)-V1)+p1(i)*n1g(i)*dnu1gdT(i)+...
        n1g(i)*RUniv;
    
    B1=[-ndot12(i); -ndot12(i)*h1f(i)*gammaN2O; 0];
    
    derivT1=A1\B1;
    
%     % Check for choked flow in vent line
%     if pAmb/p2(i)<checkSonic
%         % Flow is choked
%         mdot23(i)=A23*p2(i)*sqrt(gammaN2OKT/(RN2O*T2(i)))*...
%             (2/(gammaN2OKT+1))^((gammaN2OKT+1)/(2*(gammaN2OKT-1)));
%                 % [slug/s]
%         regime(i)=1;
%     else
%         % Flow is subsonic
%         mdot23(i)=A23*(2*gammaN2OKT/(gammaN2OKT-1)*p2(i)/nu2(i)*...
%             (1-(pAmb/p2(i))^((gammaN2OKT-1)/gammaN2OKT))*...
%             (pAmb/p2(i))^(2/gammaN2OKT))^(1/2); % [slug/s]
%         regime(i)=0;
%     end
%     
%     
%     
%     % New vent line condition
%     if p1(i)>p2(i)
%         mdot12(i)=Cd12*A12*sqrt(2*(p1(i)-p2(i))*rhof(T1(i))); % [slug/s]
%     else
%         mdot12(i)=0; % [slug/s]
%     end
%     
%     % Saturation check
%     nuSat(i)=RN2O*T2(i)/p(T2(i)); % [ft^3/slug]
%     if abs(nu2(i)-nuSat(i))>nu2(i)*.01
%         checkSat=false;
%         x2(i)=1;
%             % A quality of one for an unsaturated gas is not
%             % thermodynamically correct. However, this will be defined as
%             % one until it becomes saturated in the interest of keeping the
%             % length and index of x2 consistent with the other parameters.
%     else
%         checkSat=true
%         x2(i)=(nu2(i)-1/rhof(T2(i)))/(1/rhog(T2(i))-1/rhof(T2(i)));
%     end
    
    % Update values    
    dn2fdt(i)=derivT2(1); % [deg R/s]
    dn2gdt(i)=derivT2(2); % [mol/s]
    dT2dt(i)=derivT2(3); % [mol/s]
  
    T2(i+1)=T2(i)+dT2dt(i)*tStep; % [deg R]
    n2g(i+1)=n2g(i)+dn2gdt(i)*tStep; % [mol]
    n2f(i+1)=n2f(i)+dn2fdt(i)*tStep; % [mol]
    
    dn1fdt(i)=derivT1(1); % [deg R/s]
    dn1gdt(i)=derivT1(2); % [mol/s]
    dT1dt(i)=derivT1(3); % [mol/s]
  
    T1(i+1)=T1(i)+dT1dt(i)*tStep; % [deg R]
    n1g(i+1)=n1g(i)+dn1gdt(i)*tStep; % [mol]
    n1f(i+1)=n1f(i)+dn1fdt(i)*tStep; % [mol]
    
    mdot12(i+1)=A12*Cd12/nu1f(i)*sqrt(2*(p1(i)-p2(i))*nu1f(i)); % [slug/s] ***PLACEHOLDER***
    mdot23(i+1)=mdot12(i)*0.05; % [slug/s] ***PLACEHOLDER***
    
    % Step forward
    i=i+1;
    t(i)=t(i-1)+tStep; % [s]
    
    % Transfer mass    
    mTransfer12(i)=mdot12(i)*tStep; % [slug]
    mTransfer23(i)=mdot23(i)*tStep; % [slug]
    m1tot(i)=m1tot(i-1)-mTransfer12(i); % [slug]
    if i==1
        m2tot(i)=m2tot+mTransfer12(i)-mTransfer23(i); % [slug]
    else
        m2tot(i)=m2tot(i-1)+mTransfer12(i)-mTransfer23(i); % [slug]
    end
    
    ndot12(i)=mdot12(i)/molecMassN2O; % [mol/s]
    ndot23(i)=mdot23(i)/molecMassN2O; % [mol/s]
    
    nTransfer12(i)=mTransfer12(i)/molecMassN2O; % [mol]
    nTransfer23(i)=mTransfer23(i)/molecMassN2O; % [mol]
    
    % Compute run tank properties
    iInterp2=find(T>T2(i),1); 
    interp2=(T(iInterp2)-T2(i))/(T(iInterp2)-T(iInterp2-1)); 
    
    cp2f(i)=(cpf(iInterp2)-cpf(iInterp2-1))*interp2+cpf(iInterp2-1);
    cp2g(i)=(cpg(iInterp2)-cpg(iInterp2-1))*interp2+cpg(iInterp2-1);
    
    nu2f(i)=(nuf(iInterp2)-nuf(iInterp2-1))*interp2+cpf(iInterp2-1);
    nu2g(i)=(nug(iInterp2)-nug(iInterp2-1))*interp2+nug(iInterp2-1);
    nu2(i)=V2/m2tot(i)*molecMassN2O;
    
    p2(i)=(p(iInterp2)-p(iInterp2-1))*interp2+p(iInterp2-1);

    h2g(i)=(hg(iInterp2)-hg(iInterp2-1))*interp2+hg(iInterp2-1); % ***PLACEHOLDER***
    
    iInterp2d1=find(T2(i)-TStep<T,1); 
    interp2d1=(T(iInterp2d1)-(T2(i)-TStep))/...
        (T(iInterp2d1)-T(iInterp2d1-1));
    iInterp2d2=find(T2(i)+TStep<T,1); 
    interp2d2=(T(iInterp2d2)-(T2(i)+TStep))/...
        (T(iInterp2d2)-T(iInterp2d2-1));
    
    p2d1=(p(iInterp2d1)-p(iInterp2d1-1))*interp2d1+p(iInterp2d1-1);
    p2d2=(p(iInterp2d2)-p(iInterp2d2-1))*interp2d2+p(iInterp2d2-1);
    
    cp2fd1=(cpf(iInterp2d1)-cpf(iInterp2d1-1))*interp2d1+cpf(iInterp2d1-1);
    cp2fd2=(cpf(iInterp2d2)-cpf(iInterp2d2-1))*interp2d2+cpf(iInterp2d2-1);
    cp2gd1=(cpg(iInterp2d1)-cpg(iInterp2d1-1))*interp2d1+cpg(iInterp2d1-1);
    cp2gd2=(cpg(iInterp2d2)-cpg(iInterp2d2-1))*interp2d2+cpg(iInterp2d2-1);
    
    nu2gd1=(nug(iInterp2d1)-nug(iInterp2d1-1))*interp2d1+nug(iInterp2d1-1);
    nu2gd2=(nug(iInterp2d2)-nug(iInterp2d2-1))*interp2d2+nug(iInterp2d2-1);
    
    dp2dT(i)=(p2d2-p2d1)/(2*TStep);
    
    dcp2fdT(i)=(cp2fd2-cp2fd1)/(2*TStep);
    dcp2gdT(i)=(cp2gd2-cp2gd1)/(2*TStep);
    
    dnu2gdT(i)=(nu2gd2-nu2gd1)/(2*TStep);
    
    x2(i)=(nu2(i)-nu2f(i))/(nu2g(i)-nu2f(i));
    
    % Compute fill tank properties
    iInterp1=find(T>T1(i),1); % First index past current temperature
    interp1=(T(iInterp1)-T1(i))/(T(iInterp1)-T(iInterp1-1)); % Linear 
        % interpolation

    p1(i)=(p(iInterp1)-p(iInterp1-1))*interp1+p(iInterp1-1); % [psf]

    nu1f(i)=(nuf(iInterp1)-nuf(iInterp1-1))*interp1+nuf(iInterp1-1); % [ft^3/mol]
    nu1g(i)=(nug(iInterp1)-nug(iInterp1-1))*interp1+nug(iInterp1-1); % [ft^3/mol]
    nu1(i)=V1/m1tot(i)*molecMassN2O; % Specific volume [ft^3/slug]
    x1(i)=(nu1(i)-nu1f(i))/(nu1g(i)-nu1f(i)); % Quality in fill tank

    n1g(i)=x1(i)*n1tot; % [mol] Moles of gas in fill tank
    n1f(i)=n1tot-n1g(i); % [mol]  Moles of liquid in fill tank

    h1f(i)=(hf(iInterp1)-hf(iInterp1-1))*interp1+hf(iInterp1-1); % ***PLACEHOLDER***

    cp1f(i)=(cpf(iInterp1)-cpf(iInterp1-1))*interp1+cpf(iInterp1-1);
    cp1g(i)=(cpg(iInterp1)-cpg(iInterp1-1))*interp1+cpg(iInterp1-1);

    iInterp1d1=find(T1(i)-TStep<T,1); 
    interp1d1=(T(iInterp1d1)-(T1(i)-TStep))/...
        (T(iInterp1d1)-T(iInterp1d1-1));
    iInterp1d2=find(T1(i)+TStep<T,1); 
    interp1d2=(T(iInterp1d2)-(T1(i)+TStep))/...
        (T(iInterp1d2)-T(iInterp1d2-1));

    p1d1=(p(iInterp1d1)-p(iInterp1d1-1))*interp1d1+p(iInterp1d1-1);
    p1d2=(p(iInterp1d2)-p(iInterp1d2-1))*interp1d2+p(iInterp1d2-1);

    cp1fd1=(cpf(iInterp1d1)-cpf(iInterp1d1-1))*interp1d1+cpf(iInterp1d1-1);
    cp1fd2=(cpf(iInterp1d2)-cpf(iInterp1d2-1))*interp1d2+cpf(iInterp1d2-1);
    cp1gd1=(cpg(iInterp1d1)-cpg(iInterp1d1-1))*interp1d1+cpg(iInterp1d1-1);
    cp1gd2=(cpg(iInterp1d2)-cpg(iInterp1d2-1))*interp1d2+cpg(iInterp1d2-1);

    nu1gd1=(nug(iInterp1d1)-nug(iInterp1d1-1))*interp1d1+nug(iInterp1d1-1);
    nu1gd2=(nug(iInterp1d2)-nug(iInterp1d2-1))*interp1d2+nug(iInterp1d2-1);

    dp1dT(i)=(p1d2-p1d1)/(2*TStep);

    dcp1fdT(i)=(cp1fd2-cp1fd1)/(2*TStep);
    dcp1gdT(i)=(cp1gd2-cp1gd1)/(2*TStep);

    dnu1gdT(i)=(nu1gd2-nu1gd1)/(2*TStep);
end


% Sanity Plots
figure(1)
subplot(4,2,1)
plot(t,T2(1:i))
title('Run Tank Temperature')
subplot(4,2,3)
plot(t,p2,'r')
title('Run Tank Pressure')
subplot(4,2,5)
plot(t,nu2,'g')
title('Run Tank Specific Volume')
subplot(4,2,7)
plot(t,x2,'k')
title('Run Tank Quality')

subplot(4,2,2)
plot(t,T1(1:i))
title('Fill Tank Temperature')
subplot(4,2,4)
plot(t,p1,'r')
title('Fill Tank Pressure')
subplot(4,2,6)
plot(t,nu1,'g')
title('Fill Tank Specific Volume')
subplot(4,2,8)
plot(t,x1,'k')
title('Fill Tank Quality')

figure(2)
subplot(1,2,1)
plot(t,m2tot)
hold on
plot(t,m2tot.*x2,'g')
plot(t,m2tot-m2tot.*x2,'r')
hold off
legend('Total Mass','Mass of Gas','Mass of Liquid','Location',...
    'SouthOutside')
title('Mass of Nitrous Oxide in Run Tank')
xlabel('Time [s]')
ylabel('Mass [slug]')

subplot(1,2,2)
plot(t,m1tot)
hold on
plot(t,m1tot.*x1,'g')
plot(t,m1tot-m1tot.*x1,'r')
hold off
legend('Total Mass','Mass of Gas','Mass of Liquid','Location',...
    'SouthOutside')
title('Mass of Nitrous Oxide in Fill Tank')
xlabel('Time [s]')
ylabel('Mass [slug]')

figure(3)
plot(t,mdot12)
hold on
plot(t,mdot23,'r')
hold off
legend('Flow from Fill to Run Tank','Flow from Run Tank to Atmosphere')
title('Mass Flow Rate')