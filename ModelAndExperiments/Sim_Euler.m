function [t, y, InfCells, AgeI] = Sim_Euler(p)
%SIM_EULER solves the equation system in the function file specified by 
%   the handle p.Ex.Model at the time points p.Ex.tspan using the step size
%   p.Ex.h and the intial condition p.Ex.y0
%
%   last revised: 2018/11/21

NumSteps = length(p.Ex.tspan)-1;
InfCells = zeros(length(p.Ex.tspan), 1);
y        = [p.Ex.y0; zeros(NumSteps,length(p.Ex.y0))];

% introduce infected cell age matrix (time point x cell infection age)
AgeI     = zeros(length(p.Ex.tspan),length(p.Ex.tspan));

for Iter = 1:NumSteps
    [IntRrelI, IntRrelI_Par, IntI, IntApoI, AgeI(Iter,:)] = Sim_Integrals(y(1:Iter,:), p);
    y(Iter+1,:) = max(zeros(size(p.Ex.y0)), ...
        y(Iter,:) + p.Ex.h*feval(p.Ex.ModelName, p.Ex.tspan(Iter), y(Iter,:),...
                                 IntRrelI, IntRrelI_Par, IntI, IntApoI, p));
    y(Iter+1, y(Iter+1,:)<1e-20) = 0;
    InfCells(Iter) = IntI;
    
    %prevent overloading of high-affinity binding sites due to numerical inaccuracies
    if(p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) - y(Iter+1,p.Idx.VAttHi) < 0)

        %save attached virus concentration
        a = y(Iter+1, p.Idx.VAttHi);
        
        %calculate attached virus concentration that saturates all binding sites
        y(Iter+1, p.Idx.VAttHi) = p.Ex.BtotHi*sum(y(Iter+1,p.Idx.T|p.Idx.Ta));

        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.V)      = y(Iter+1, p.Idx.V) + a-y(Iter+1, p.Idx.VAttHi);
    end
    
    %prevent overloading of low-affinity binding sites due to numerical inaccuracies
    if(p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta)) - y(Iter+1,p.Idx.VAttLo) < 0)

        %save attached virus concentration
        a = y(Iter+1, p.Idx.VAttLo);
        
        %calculate attached virus concentration that saturates all binding sites
        y(Iter+1, p.Idx.VAttLo) = p.Ex.BtotLo*sum(y(Iter+1,p.Idx.T|p.Idx.Ta));
        
        %transfer the excessive virions back to the non-attached state
        y(Iter+1, p.Idx.V)      = y(Iter+1, p.Idx.V) + a-y(Iter+1, p.Idx.VAttLo);
    end
    
end

[~, ~, InfCells(Iter+1), ~] = Sim_Integrals(y, p);
t = p.Ex.tspan;
