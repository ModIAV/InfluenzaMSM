function result = Opt_Optimization
%OPT_OPTIMIZATION provides the interface between the optimization solver and
%   the objective function
%
%   last revised: 2018/11/21

global p
global d

%% Optimization solver options
    %% Fminsearch and Fmincon
    p.FminOptions = optimset('MaxIter',     p.MaxIter,...
                             'MaxFunEvals', p.MaxEval,...
                             'Display',     'iter');
                       
    %% fSSm
    p.fSSmOptions.log_var      = 1:p.In.nOpt+p.Ex.nOpt;
    p.fSSmOptions.maxtime      = p.MaxTime;            % maximal time limit in seconds (e.g. 1 day = 60*60*24)
    p.fSSmOptions.maxeval      = p.MaxEval;            % maximal number of function evaluations (e.g. 1000000)
    p.fSSmOptions.local.solver = 0;
    
    %% SRES
    p.SRESOptions.lambda = 100;                        % population size (number of offspring) (100 to 200)
    p.SRESOptions.G      = 50;                         % maximum number of generations
    p.SRESOptions.mu     = 1/7*p.SRESOptions.lambda;   % parent number (mu/lambda usually 1/7)
    p.SRESOptions.pf     = 0.45;                       % pressure on fitness in [0 0.5] try around 0.45
    p.SRESOptions.varphi = 1;                          % expected rate of convergence (usually 1)
    
result = [];    
    
%% Prepare optimization of the intracellular model
% Find position of initial conditions subject to optimization in p.In.Ic
if(isempty(p.In.IcOpt))
    p.In.PosIcOpt = [];
else
    p.In.PosIcOpt = zeros(length(p.In.IcOpt),1);
    for CountIc=1:length(p.In.IcOpt)
        if(sum(strcmp(p.In.IcOpt(CountIc), p.In.StateNames))==0)
            disp(char(strcat('Warning(Optimization): "', p.In.IcOpt(CountIc), '" is no initial condition in the intracellular model.')));
            return
        else
            p.In.PosIcOpt(CountIc) = find(strcmp(p.In.IcOpt(CountIc), p.In.StateNames));
        end
    end
end

% Set optimization bounds and initial guess for the intracellular model
p.In.x_L = zeros(p.In.nOpt,1);
p.In.x_U = zeros(p.In.nOpt,1);
p.In.x0  = zeros(p.In.nOpt,1);
for CountOptPa=1:length(p.In.PaOpt)
    if(sum(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))==0)
        disp(strcat('Warning(Optimization): "', p.In.PaOpt(CountOptPa),...
            '" is no parameter in the intracellular model.'));
        return
    end
    p.In.x_L(CountOptPa) = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))/p.In.PaDev;
    p.In.x_U(CountOptPa) = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames))*p.In.PaDev;
    p.In.x0(CountOptPa)  = p.In.ParaValues(strcmp(p.In.PaOpt(CountOptPa), p.In.ParaNames));
end
if(~isempty(p.In.IcOpt))
    for CountOptIc=1:length(p.In.IcOpt)
        p.In.x_L(length(p.In.PaOpt)+CountOptIc) = p.In.Ic(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames))/p.In.IcDev;
        p.In.x_U(length(p.In.PaOpt)+CountOptIc) = p.In.Ic(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames))*p.In.IcDev;
        p.In.x0(length(p.In.PaOpt)+CountOptIc)  = p.In.Ic(strcmp(p.In.IcOpt(CountOptIc), p.In.StateNames));
    end
end

% Find position of states subject to optimization in the intracellular model
for CountFit=1:length(p.In.States2Fit)
    if(sum(strcmp(p.In.States2Fit(CountFit), [p.In.StateNames; p.In.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.In.States2Fit(CountFit), '" is no state or variable in the intracellular model.'));
        return
    elseif(sum(strcmp(p.In.States2Fit(CountFit), d.In.ExpStateNames))==0)
        disp(strcat('Warning(Optimization): "', p.In.States2Fit(CountFit), '" has not been measured.'));
        return
    else
        p.In.PosState2Fit(CountFit, 1) = find(strcmp(p.In.States2Fit(CountFit), [p.In.StateNames; p.In.VariableNames]));
        p.In.PosState2Fit(CountFit, 2) = find(strcmp(p.In.States2Fit(CountFit), d.In.ExpStateNames));
    end
end

%% Prepare optimization of the intracellular model for coupling (version w/o virus entry)
% Find position of initial conditions subject to optimization in p.InToEx.Ic
p.InToEx.PosIcOpt = p.In.PosIcOpt - (length(p.In.Ic)-length(p.InToEx.Ic));    % positions of optimized intracellular IC using diff(# of states) between In and InToEx
p.InToEx.PosIcOpt(p.InToEx.PosIcOpt==1) = 0;
p.InToEx.PosIcOpt(p.InToEx.PosIcOpt==-4) = 1;                                 % if Vex(0) is changed, adjust VpCyt(0) accordingly, disregard IC not present in InToEx
p.InToEx.PosIcOpt(p.InToEx.PosIcOpt<=0) = 0;

%% Prepare optimization of the extracellular model
% Find position of initial conditions subject to optimization in p.Ex.Ic
if(isempty(p.Ex.IcOpt)) 
    p.Ex.PosIcOpt = [];
else
    p.Ex.PosIcOpt = zeros(length(p.Ex.IcOpt),1);
    for CountIc=1:length(p.Ex.IcOpt)
        if(sum(strcmp(p.Ex.IcOpt(CountIc), p.Ex.StateNames))==0)
            disp(char(strcat('Warning(Optimization): "', p.Ex.IcOpt(CountIc), '" is no initial condition in the extracellular model.')));
            return
        else
            p.Ex.PosIcOpt(CountIc) = find(strcmp(p.Ex.IcOpt(CountIc), p.Ex.StateNames));
        end
    end
end

% Set optimization bounds and initial guess for the extracellular model
p.Ex.x_L = zeros(p.Ex.nOpt,1);
p.Ex.x_U = zeros(p.Ex.nOpt,1);
p.Ex.x0  = zeros(p.Ex.nOpt,1);
for CountOptPa=1:length(p.Ex.PaOpt)
    if(sum(strcmp(p.Ex.PaOpt(CountOptPa), p.Ex.ParaNames))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.PaOpt(CountOptPa),...
            '" is no parameter in the extracellular model.'));
        return
    end
    p.Ex.x_L(CountOptPa) = p.Ex.(char(p.Ex.PaOpt(CountOptPa)))/p.Ex.PaDev;
    p.Ex.x_U(CountOptPa) = p.Ex.(char(p.Ex.PaOpt(CountOptPa)))*p.Ex.PaDev;
    p.Ex.x0(CountOptPa)  = p.Ex.(char(p.Ex.PaOpt(CountOptPa)));
end
if(~isempty(p.Ex.IcOpt))
    for CountOptIc=1:length(p.Ex.IcOpt)
        p.Ex.x_L(length(p.Ex.PaOpt)+CountOptIc) = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames))/p.Ex.IcDev;
        p.Ex.x_U(length(p.Ex.PaOpt)+CountOptIc) = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames))*p.Ex.IcDev;
        p.Ex.x0(length(p.Ex.PaOpt)+CountOptIc)  = p.Ex.y0(strcmp(p.Ex.IcOpt(CountOptIc), p.Ex.StateNames));
    end
end


% Find position of states subject to optimization in the extracellular model
for CountFit=1:length(p.Ex.States2Fit)
    if(sum(strcmp(p.Ex.States2Fit(CountFit), [p.Ex.StateNames p.Ex.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.States2Fit(CountFit), '" is no state or variable in the extracellular model.'));
        return
    elseif(sum(strcmp(p.Ex.States2Fit(CountFit), d.Ex.ExpStateNames))==0)
        disp(strcat('Warning(Optimization): "', p.Ex.States2Fit(CountFit), '" has not been measured.'));
        return
    else
        p.Ex.PosState2Fit(CountFit, 1) = find(strcmp(p.Ex.States2Fit(CountFit), [p.Ex.StateNames p.Ex.VariableNames]));
        p.Ex.PosState2Fit(CountFit, 2) = find(strcmp(p.Ex.States2Fit(CountFit), d.Ex.ExpStateNames));
    end
end

%% Define weighting matrix for intracellular model
p.In.W = 1./repmat(max(d.In.ExpStateValues(:, p.In.PosState2Fit(:,2))), length(d.In.Time), 1); %normalize to maximum of measurement
p.In.W(isnan(d.In.ExpStateValues(:, p.In.PosState2Fit(:,2)))) = 0;

%% Define weighting matrix for extracellular model
p.Ex.W = 1./repmat(max(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2))), length(d.Ex.Time), 1); %normalize to maximum of measurement
p.Ex.W(isnan(d.Ex.ExpStateValues(:, p.Ex.PosState2Fit(:,2)))) = 0;

%% Prepare vector to realize offset in RNA concentrations for intracellular model
% as free viral RNAs in the seed virus might adhere to cells and thus cause
% a constant offset in measurements we implement such an offset in the
% objective function (*_ObjFun.m)
RnaNames = {'Rm5', 'RcSegment', 'RvSegment'};
p.In.OffsetIndex = zeros(length(RnaNames), 2);
for CountRnas = 1:length(RnaNames)
    p.In.OffsetIndex(CountRnas, 1) = find(strcmp(RnaNames(CountRnas), [p.In.StateNames; p.In.VariableNames]));
    p.In.OffsetIndex(CountRnas, 2) = find(strcmp(RnaNames(CountRnas), d.In.ExpStateNames));
end

%% Check whether parameter occuring in both models are subject to optimization
p.In.IdxKfus  = false(1,p.In.nPaOpt);
p.In.IdxKfus  = strcmp('kFus', p.In.PaOpt);

if(sum(strcmp('kFus', p.Ex.PaOpt)))
    disp('Warning(Optimization): If kFus should be optimized include it only in p.In.PaOpt instead of p.Ex.PaOpt.');
    return
end

ParasInBothModels = {'NHi', 'NLo', 'kAtt', 'kEqHi', 'kEqLo', 'kEn'}; 
for CountOptPa=1:length(p.In.PaOpt)
    if(sum(strncmp(p.In.PaOpt(CountOptPa), ParasInBothModels, 3)))
        disp('Warning(Optimization): Optimization of entry parameters which occur in both the intracellular and extracellular model is not implemented (except for kFus).');
        return
    end
end
for CountOptPa=1:length(p.Ex.PaOpt)
    if(sum(strncmp(p.Ex.PaOpt(CountOptPa), ParasInBothModels, 3)))
        disp('Warning(Optimization): Optimization of entry parameters which occur in both the intracellular and extracellular model is not implemented (except for kFus).');
        return
    end
end

%% Start parameter estimation

switch lower(p.Solver)
    case {'fminsearch'}
        %% Fminsearch
        [xbest, result.ObjFun] = ...
            fminsearch(@Opt_ObjFun, [p.In.x0; p.Ex.x0], p.FminOptions);
        result.In.xbest = xbest(1:p.In.nOpt);
        result.Ex.xbest = xbest(p.In.nOpt+1:end);

    case {'fmincon'}
        [xbest, result.ObjFun] = ...
            fmincon(@Opt_ObjFun, [p.In.x0; p.Ex.x0],[],[],[],[],...
            [p.In.x_L; p.Ex.x_L],[p.In.x_U; p.Ex.x_U],[],p.FminOptions);
        result.In.xbest = xbest(1:p.In.nOpt);
        result.Ex.xbest = xbest(p.In.nOpt+1:end); 

    case {'fssm'}
        %% fSSm
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

        p.In.x_U(strcmp('FPar',[p.In.PaOpt p.In.IcOpt])) = 1; % set upper threshold for FPar <= 1
        
        p.problem.f       = 'Opt_ObjFun';
        p.problem.x_0     = [p.In.x0;  p.Ex.x0];
        p.problem.x_L     = [p.In.x_L; p.Ex.x_L];
        p.problem.x_U     = [p.In.x_U; p.Ex.x_U];
        p.problem.int_var = 0;
        p.problem.bin_var = 0;
        p.problem.c_L     = [];
        p.problem.c_U     = [];
        p.problem.N       = numel(d.In.Time) + numel(d.Ex.Time);

        [result.fSSm] = fssm_kernel(p.problem, p.fSSmOptions);

        xbest  = result.fSSm.xbest;

        result.In.xbest  = xbest(1:p.In.nOpt);
        result.Ex.xbest  = xbest(p.In.nOpt+1:end);
        result.ObjFunVal = result.fSSm.fbest;

    otherwise
        warning(char(['Optimization - No optimization algorithm "',p.solver,'" found!']));
        result = [];
end

