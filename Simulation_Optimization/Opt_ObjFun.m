function ObjFunVal = Opt_ObjFun(x)
%OPT_OBJFUN determines the deviation of measurements and simulation by using W 
%   as the weighting matrix and returns a scalar for optimization routines
%
%   last revised: 2018/11/27

global p
global d
clear mex

xIn = x(1:p.In.nOpt);
xEx = x(p.In.nOpt+1:end);

%% Intracellular model
%adopt new initial conditions
p.In.Ic(p.In.PosIcOpt) = xIn(p.In.nPaOpt+1:p.In.nOpt);
%simulate the model with new parameters and initial conditions
try
    resIn = SBPDsimulate('InModel_MexFile', d.In.Time, p.In.Ic, p.In.PaOpt, xIn(1:p.In.nPaOpt));   
catch ModelCompileError
    InModel_tmp = SBparameters(p.InModel, p.In.PaOpt, xIn(1:p.In.nPaOpt));
    resIn       = SBsimulate(InModel_tmp, 'ode23s', d.In.Time, p.In.Ic); 
end
SimulationIn = [resIn.statevalues, resIn.variablevalues];

% To account for free viral RNAs in the seed virus which might adhere to
% cells and cause a constant offset in measurements, we add the 0 hpi
% simulation value of each RNA species to the subsequent time points
SimulationIn(:,p.In.OffsetIndex(:,1)) = SimulationIn(:,p.In.OffsetIndex(:,1)) +...
    repmat(d.In.ExpStateValues(1,p.In.OffsetIndex(:,2)), length(resIn.time), 1);

SimValuesIn      =        SimulationIn(:,p.In.PosState2Fit(:,1));
MeasuredValuesIn = d.In.ExpStateValues(:,p.In.PosState2Fit(:,2));

% To yield an HA titer cells need to produce more than 2e7 virions/ml,
% which means that with a cell concentration of roughly 1.2e6/ml each cell
% must at least produce 16.7 virions to yield a measureable HA titer
a = SimValuesIn(:,strcmp('Vrel', p.In.States2Fit));
a(a<16.7) = 0;
SimValuesIn(:,strcmp('Vrel', p.In.States2Fit)) = a;

DevStatesIn = SimValuesIn-MeasuredValuesIn;
DevStatesIn(isnan(MeasuredValuesIn)) = 0;
ObjFunValIn = sum(sum((p.In.W .* DevStatesIn).^2).*p.In.RelStateWeight);

%% Simulate the intracellular model for coupling (version w/o virus entry)
% adopt new initial conditions
if ( p.InToEx.PosIcOpt(1) == 1 )
    p.InToEx.Ic(1) = max([8 8*xIn(p.In.nPaOpt+1)*p.In.ParaValues(strcmp('Ffus', p.In.ParaNames))]);
    if ( length(p.InToEx.PosIcOpt) > 2 )
        p.InToEx.Ic(p.InToEx.PosIcOpt(2:end)) = xIn(p.In.nPaOpt+2:p.In.nOpt);
    end
else        
    p.InToEx.Ic(p.InToEx.PosIcOpt) = xIn(p.In.nPaOpt+1:p.In.nOpt);         
end

% get new parameters        
xInToEx = xIn(1:p.In.nPaOpt);
try
    resInToEx = SBPDsimulate('InModelToEx_MexFile', p.Ex.tspan, p.InToEx.Ic,...
        p.In.PaOpt(~p.In.IdxKfus),xInToEx(~p.In.IdxKfus));
catch ModelCompileError
    InToExModel_tmp = SBparameters(p.InToExModel, p.In.PaOpt(~p.In.IdxKfus),xInToEx(~p.In.IdxKfus));
    resInToEx       = SBsimulate(InToExModel_tmp, 'ode23s', p.Ex.tspan, p.InToEx.Ic); 
end

p.InToEx.Rrel     = resInToEx.reactionvalues(:,strcmp('rRel', resInToEx.reactions));
p.InToEx.Rrel_Par = resInToEx.reactionvalues(:,strcmp('rRel_par', resInToEx.reactions));

%% Extracellular model
%adopt new initial conditions
p.Ex.y0(p.Ex.PosIcOpt) = xEx(p.Ex.nPaOpt+1:p.Ex.nOpt);

%adopt parameter changes
for CountPara = 1:p.Ex.nPaOpt
    p.Ex.(char(p.Ex.PaOpt(CountPara))) = xEx(CountPara);
end
if(sum(strcmp('kFus', p.In.PaOpt)));   p.Ex.kFus  = xIn(p.In.IdxKfus);  end

%apply parameter constraint
if(p.Ex.Finf>1 || p.Ex.Finf<0)
    ObjFunVal = 1e8;
    return
end

%modify apoptosis rates and integrals
[p.Ex.ApoI, p.Ex.IntKapoi] = getApoI(p);

p.Ex.RrelIntKapoi     = p.InToEx.Rrel'.*p.Ex.IntKapoi;
p.Ex.RrelIntKapoi_Par = p.InToEx.Rrel_Par'.*p.Ex.IntKapoi;

%simulate the system with new parameters and initial conditions
[resEx.time, resEx.statevalues, resEx.variablevalues(:, p.Idx.I)] = Sim_Euler(p);

resEx.variablevalues(:,strcmp('Log10Tcid50', p.Ex.VariableNames)) =...
    log10(max(1, resEx.statevalues(:,p.Idx.V)));
resEx.variablevalues(:,strcmp('Log10Ha', p.Ex.VariableNames)) =...
    log10(max(1, resEx.statevalues(:,p.Idx.Par))); 

%calculate the deviation between simulation and measurements
SimulationEx     = [resEx.statevalues resEx.variablevalues];
SimValuesEx      = SimulationEx(p.Ex.TimeIndexData, p.Ex.PosState2Fit(:,1));
MeasuredValuesEx =           d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2));

DevStatesEx = SimValuesEx-MeasuredValuesEx;
DevStatesEx(isnan(MeasuredValuesEx)) = 0;
ObjFunValEx = sum(sum((p.Ex.W .* DevStatesEx).^2).*p.Ex.RelStateWeight);

%% Combine both objective function values
ObjFunVal = ObjFunValIn/(numel(MeasuredValuesIn)-sum(sum(isnan(MeasuredValuesIn))))...
           +ObjFunValEx/(numel(MeasuredValuesEx)-sum(sum(isnan(MeasuredValuesEx))));

% ObjFunVal = ObjFunVal*1e6;
