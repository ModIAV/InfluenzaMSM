function [IntRrelI, IntRrelI_Par, IntI, IntApoI, AgeI] = Sim_Integrals(y, p)
%SIM_INTEGRALS calculates the integrals of virus production, infected
%   cells and the apoptosis rate of infected cells over all cell ages tau.
%
%   last revised: 2018/11/21

n = size(y,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define integrals for tau=t, cells that were (already) infected at t=0
IntI         = p.Ex.I0*p.Ex.IntKapoi(n);
IntRrelI     = IntI*p.Ex.RrelIntKapoi(n);
IntRrelI_Par = IntI*p.Ex.RrelIntKapoi_Par(n); %
IntApoI      = IntI*p.Ex.ApoI(n);

% infected cell age vector, collects how many cells have been infected for 
% how long at the current time point
AgeI = zeros(1,length(p.Ex.tspan));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InfForce = zeros(n, 1);

TargetCells = y(:,p.Idx.T)+y(:,p.Idx.Ta);
InfForce(TargetCells>0) = p.Ex.Finf*p.Ex.kFus*y(TargetCells>0,p.Idx.VEn)...
    .*y(TargetCells>0,p.Idx.T)./TargetCells(TargetCells>0);

for j=1:n   
    IntI         = IntI         + InfForce(n+1-j)        *p.Ex.IntKapoi(j) * p.Ex.h;  
    IntRrelI     = IntRrelI     + InfForce(n+1-j)    *p.Ex.RrelIntKapoi(j) * p.Ex.h;
    IntRrelI_Par = IntRrelI_Par + InfForce(n+1-j)*p.Ex.RrelIntKapoi_Par(j) * p.Ex.h;
    IntApoI      = IntApoI      + InfForce(n+1-j)        *p.Ex.IntKapoi(j) * p.Ex.h * p.Ex.ApoI(j);
end

