function [MOI,MOIsd] = calculateMOI(V0,H0,V0sd,H0sd)

MOI = V0./H0;
MOIsd = sqrt(sum([(-V0./(H0.^2)).*H0sd (1./H0).*V0sd].^2,2));