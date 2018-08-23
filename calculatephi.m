function [phi,phisd] = calculatephi(V0,V1,H0,H1,V0sd,V1sd,H0sd,H1sd)
%Length of time before washout (hours)
T = .25;
%vector of fraction of hosts lost
fh = (H0-H1)./H0;
%Standard deviation with propagation of uncertainty
fhsd = sqrt(sum([(H1./(H0.^2)).*H0sd (-1./H0).*H1sd].^2,2));
%mean loss
fhbar = mean(fh);
%propagate uncertainty -- treat as mixture distribution for averaged
%quantity
fhbarsd = std_biotech(fh,fhsd);
%Need stddev of probability of infection, this is done by Taylor expansion
%common denominator of Taylor expansion
commondenom = (1./(1-fhbar));
Pinfsd = sqrt(sum([commondenom.*(-V1./(V0.^2)).*V0sd commondenom.*(1./V0).*V1sd ...
    commondenom.^2.*(1+(V1-V0)./V0).*fhbarsd].^2,2));
Pinf = commondenom.*(1+(V1-V0)./V0);

%adsorption rate function
phi = -log(1-Pinf)./T./H0;
%propagate uncertainty
phisd = sqrt(sum([(1./((1-Pinf).*H0.*T)).*Pinfsd log(1-Pinf)./(T.*(H0).^2).*H0sd].^2,2));