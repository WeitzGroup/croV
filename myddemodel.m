function [sol] = myddemodel(lattime,Htot,Hinf,numboxcar)
%Simulated dynamics where number of infected host subcompartments (boxcars) is variable. 
    
%initial conditions
V0 = 0;
H0 = Htot-Hinf;
I10 = Hinf;


initialvec = [V0; H0; 0; I10; zeros(numboxcar-1,1)];
options = ddeset('InitialY',initialvec);
%delay differential equation
numstate = numboxcar+3;
sol = dde23(@ddex1de,lattime,zeros(numstate,1),[0, 24-lattime],options);

function dydt = ddex1de(t,y,Z)
%Generate delay differential equation model

%define variables
vlag = Z(1);
hlag = Z(2);
V = y(1);
H = y(2);
I0 = y(3);
I1 = y(4);
Iend = y(end);
%set parameters
phi = 1*10^(-6);
%phipar = 0;
r = 0.2;
K = 5.41*10^6;
burstsize = 450;
%scale transitions through subcompartments
lysisrate = numboxcar/3;

%set functions for lagged and non-lagged expressions
phifun = @(currv,currh) phi*currv*currh;
Ifun = @(currv,currh) phi.*currh.*currv.*phifun(currv,currh);
%Delay differential equations
dVdt = burstsize*lysisrate*Iend-Ifun(V,H);
dHdt = r.*(H).*(1-(H+I0)./K) - Ifun(V,H);
dI0dt = Ifun(V,H) - Ifun(vlag(1),hlag(1));
dI1dt = Ifun(vlag(1),hlag(1)) - lysisrate*I1;
dydt = [dVdt; dHdt; dI0dt; dI1dt; zeros(numboxcar-1,1)];

%loop around subcompartments and add DDE terms
boxidx = 4;
for jj = 1:numboxcar-1
    boxidx = boxidx+1;
    dydt(boxidx) = lysisrate*(y(boxidx-1)-y(boxidx));
end
end

end


