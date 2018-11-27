%%% --------------------------------------------------------------- %%%
%%% Simulation and Optimization of a Multi-scale model of influenza %%% 
%%% virus infection in cell cultures                                %%%
%%% --------------------------------------------------------------- %%%

% authors: Daniel Ruediger & Stefan Heldt
% last revised: 2018/11/27

%% Housekeeping commands I
clear all; close all; clear mex; clc;
global p d

delete InModel_MexFile.mexa64           % deletion of Mex-files prevents occasional program crashes
delete InModelToEx_MexFile.mexa64

ModelPath = pwd; Index = find(filesep==ModelPath);
ModelPath = strcat(ModelPath(1:Index(end)), 'ModelAndExperiments/');
p.SearchPath = path; path(ModelPath, path);
clear Index
    
p.Info.StartTime = clock;
p.Ex = Sim_ParameterDeclaration; % import basic extracellular parameters

%% Main options
p.Task    = 'optimize';  % tasks: 'simulate', 'optimize'
p.SimMode = '';          % modes: ''        (extended model) 
                         %        'noInhib' (extended model w/o mRNA synthesis inhibition)

p.Ex.Moi  = 73; 	% multiplicity of infection
                    % - MOI 73 simulates the model with the parameters given in the paper
                    % - MOI 3 and 1e-4 show the model prediction compared to experimental data
                    % - other MOIs simulate a model prediction for the specific conditions

p.Apo.Increase = 'logistic';      % choose apoptosis rate increase function:
                                  % 'none', 'linear', 'logistic', 'normal_distribution', 'gompertz', 'hill*' (*=1-4), 

%% Optimization options
p.Solver  = 'fminsearch';      %optimization solver: 'fSSm', 'fminsearch', 'fmincon'
p.MaxTime = 2*60*60;           %maximum run time of optimization solver in seconds (fSSm)
p.MaxIter = 1e3;               %maximum iterations
p.MaxEval = 1e6;               %maximum function evaluations

% intracellular model
p.In.PaOpt = {'kFus','kSynV', 'kSynC', 'kSynM', 'kDegM','kBindM1','kRelRed', 'kRel', 'KvRel'};     %parameters subject to optimization in the intracellular model
if ( ~strcmp('noInhib',p.SimMode) )    % add parameter Kr to the estimated paramters for model with mRNA synthesis inhibition
    p.In.PaOpt{end+1} = 'Kr';
end
p.In.IcOpt = {'FPar'};    %initial conditions subject to optimization in the intracellular model, set Vex as the first variable
p.In.PaDev = 50;                 %upper and lower bounds of parameters deviate by */ PaDev from initial parameter guess
p.In.IcDev = 10;                 %upper and lower bounds of IC deviate by */ IcDev from IC

p.In.States2Fit = {'Rm5', 'RcSegment', 'RvSegment', 'Vrel', 'Prel'}; %states subject to fitting for the intracellular model
p.In.RelStateWeight = ones(size(p.In.States2Fit));     %relative weight of the deviation of each fitted state in the objective function

% extracellular model
p.Ex.PaOpt = {'kLys', 'kDegV', 'kApoT', 'kApoI', 'tau_Apo', 'v_Apo'}; 
p.Ex.IcOpt = {};        %initial conditions subject to optimization in the extracellular model
p.Ex.PaDev = 50;                 %upper and lower bounds of parameters deviate by */ PaDev from initial parameter guess
p.Ex.IcDev = 10;               %upper and lower bounds of IC deviate by */ IcDev from IC

p.Ex.States2Fit = {'T', 'I', 'Ia', 'Ta', 'V', 'Par'};             %states subject to fitting for the extracellular model
p.Ex.RelStateWeight = ones(size(p.Ex.States2Fit));         %relative weight of the deviation of each fitted state in the objective function

p.In.nPaOpt = length(p.In.PaOpt);                          %number of parameters subject to optimization
p.In.nOpt   = length(p.In.PaOpt)+length(p.In.IcOpt);       %number of parameters and initial conditions subject to optimization 
p.Ex.nPaOpt = length(p.Ex.PaOpt);                          %number of parameters subject to optimization
p.Ex.nOpt   = length(p.Ex.PaOpt)+length(p.Ex.IcOpt);       %number of parameters and initial conditions subject to optimization 

%% Final parameters and initial conditions

if ( ~strcmp(p.SimMode, 'noInhib') )
    % Intracellular model
    p.In.ChangedParaNames  = { 'kSynV', 'kSynC', 'kSynM', 'kDegM', 'kBindM1', 'kRelRed', 'kRel', 'KvRel', 'Kr'};     %names  of parameters that ought to be changed in the intracellular model compared to the parameter declaration file
    p.In.ChangedParaValues = [   8.4,     0.8,    1.8e5,   0.33,     9e-7,      0.05,     1270,   1250,   1.1e7];     

    p.In.ChangedIcNames  = { 'FPar'};    %names  of initial conditions that ought to be changed in the intracellular model compared to the parameter declaration file
    p.In.ChangedIcValues = 3.4e-2;

    % Extracellular model
    p.Ex.ChangedParaNames  = { 'kLys', 'kDegV', 'kApoT', 'kApoI', 'v_Apo', 'tau_Apo'};      %names  of parameters that ought to be changed in the extracellular model compared to the parameter declaration file
    p.Ex.ChangedParaValues = [  9.4e-3, 1.15e-2,  7e-3,    0.11,   0.77,    19.8];

    % Parameters included on both scales
    p.ChangedParaNames     = { 'kFus'};     %names  of parameters that ought to be changed in both the intracellular and extracellular model (mainly parameters for virus entry)    
    p.ChangedParaValues    = 0.31;
else
    % Intracellular model
    p.In.ChangedParaNames  = { 'kSynV', 'kSynC', 'kSynM', 'kDegM', 'kBindM1', 'kRelRed', 'kRel', 'KvRel'};     %names  of parameters that ought to be changed in the intracellular model compared to the parameter declaration file
    p.In.ChangedParaValues = [   8.3,     0.9,    1.7e5,   0.63,    1.1e-6,     0.05,    1270,    1250];     
    
    p.In.ChangedIcNames  = { 'FPar'};    %names  of initial conditions that ought to be changed in the intracellular model compared to the parameter declaration file
    p.In.ChangedIcValues = 3.5e-2;
    
    % Extracellular model
    p.Ex.ChangedParaNames  = {'kLys', 'kDegV', 'kApoT', 'kApoI', 'v_Apo', 'tau_Apo'};      %names  of parameters that ought to be changed in the extracellular model compared to the parameter declaration file
    p.Ex.ChangedParaValues = [ 9.3e-3, 1.15e-2,  7e-3,   0.11,   0.76,    19.8];
    
    % Parameters included on both scales
    p.ChangedParaNames     = {'kFus'};     %names  of parameters that ought to be changed in both the intracellular and extracellular model (mainly parameters for virus entry)    
    p.ChangedParaValues    =   0.31;      %values of parameters that ought to be changed in both the intracellular and extracellular model (mainly parameters for virus entry)   
end

if (strcmpi(p.Apo.Increase(1:4),'hill')) %extract hill coefficient from options
    p.Ex.ChangedParaValues(5) = str2double(p.Apo.Increase(end));
end

p.In.h  = 0.02;                 %h  step-size of the intracellular model
p.Ex.h  = 0.05;                 %h  step-size used to solve the extracellular model with Euler's-Method

%% Check toolbox availability
if ( ~exist('SBmodel') )
    error('For model simulation and optimization the Systems Biology Toolbox (Schmidt et al., 2006, J Glob Opt) is required. You can download it at http://www.sbtoolbox2.org/main.php?display=download&menu=download');
end

%% Load model and apply parameter changes
if ( ~strcmp(p.SimMode, 'noInhib') )
    p.In.ModelName     = 'Sim_IntracellularLevel_FullModel.txt';      % name of intracellular model
    p.InToEx.ModelName = 'Sim_IntracellularLevel_ReducedModel.txt';   % name of intracellular model for coupling to the extracellular model
else
    p.In.ModelName     = 'Sim_IntracellularLevel_FullModel_noInhib.txt';      % name of intracellular model
    p.InToEx.ModelName = 'Sim_IntracellularLevel_ReducedModel_noInhib.txt';   % name of intracellular model for coupling to the extracellular model
end

%load intracellular models
InModel     = SBmodel(fullfile(ModelPath, p.In.ModelName));
InToExModel = SBmodel(fullfile(ModelPath, p.InToEx.ModelName));

%apply parameter changes to the intracellular model
InModel     = SBparameters(InModel, horzcat(p.ChangedParaNames, p.In.ChangedParaNames), [p.ChangedParaValues p.In.ChangedParaValues]);
InToExModel = SBparameters(InToExModel,  p.In.ChangedParaNames, p.In.ChangedParaValues);

%apply parameter changes to the extracellular model
h.Ex.ChangedParaNames  = horzcat(p.ChangedParaNames, p.Ex.ChangedParaNames);
h.Ex.ChangedParaValues = [p.ChangedParaValues, p.Ex.ChangedParaValues];
for CountPara = 1:length(h.Ex.ChangedParaNames)
    if(isfield(p.Ex, h.Ex.ChangedParaNames(CountPara)))
        p.Ex.(char(h.Ex.ChangedParaNames(CountPara))) = h.Ex.ChangedParaValues(CountPara);
    else
        disp(strcat('Warining(Main): Parameter "', h.Ex.ChangedParaNames(CountPara),...
            '" is no parameter in the model.'));
        return
    end
end

%% load experiments and merge with model
% load respective experimental data
if ( p.Ex.Moi == 73 )
    p.In.ExpName   = 'Data_MOI73_intracellular';           % experiment name for the intracellular level
    p.Ex.ExpName   = 'Data_MOI73_extracellular';           % experiment name for the extracellular level
else   
    d.In.Time = 24;
    if ( p.Ex.Moi == 3 )
        p.Ex.ExpName   = 'Data_MOI3_titers';           % experiment name for the extracellular level       
    elseif ( p.Ex.Moi == 1e-4 )
        p.Ex.ExpName   = 'Data_MOI00001_titers';       % experiment name for the extracellular level      
    else
        p.Ex.ExpName   = 'Data_MOI3_titers';           % load experiment file anyway to obtain basic simulation info
    end
end

%intracellular model
if ( p.Ex.Moi == 73 )
    InMeasurements = SBmeasurement(fullfile(ModelPath, strcat(p.In.ExpName, '.csv')));
    InExperiment   =  SBexperiment(fullfile(ModelPath, strcat(p.In.ExpName, '.exp')));
    [d.In.Time, d.In.ExpStateNames, d.In.ExpStateValues] = SBmeasurementdata(InMeasurements);
    InModel = SBmergemodexp(InModel, InExperiment);
end

%extracellular model
ExMeasurements = SBmeasurement(fullfile(ModelPath, strcat(p.Ex.ExpName, '.csv')));
ExExperiment   =  SBexperiment(fullfile(ModelPath, strcat(p.Ex.ExpName, '.exp')));
[d.Ex.Time, d.Ex.ExpStateNames, d.Ex.ExpStateValues] = SBmeasurementdata(ExMeasurements);

p.In.tspan = 0:p.In.h:max(d.In.Time);  %h  time points for the simualtion of the intracellular model
p.Ex.tspan = 0:p.Ex.h:max(d.Ex.Time);  %h  time points used to solve the extracellular model by Euler's-Method

for CountTime=1:length(d.Ex.Time)
    if(isempty(find(p.Ex.tspan==d.Ex.Time(CountTime), 1)))
        disp(strcat('Warining(Main): Choose p.Ex.h such that a vector 0:p.Ex.h:max([p.Ex.Time; d.In.Time]) contains all measurement time points.'));
        return
    else
        p.Ex.TimeIndexData(CountTime) = find(p.Ex.tspan==d.Ex.Time(CountTime), 1);
    end
end

%% Prepare simulation of the intracellular model
p.In.VariableNames                        = SBvariables(InModel);
p.InToEx.VariableNames                    = SBvariables(InToExModel);
[p.In.StateNames, ~, p.In.Ic]             =     SBstates(InModel);
[p.InToEx.StateNames, ~, p.InToEx.Ic]     =     SBstates(InToExModel);
[p.In.ParaNames, p.In.ParaValues]         = SBparameters(InModel);
[p.InToEx.ParaNames, p.InToEx.ParaValues] = SBparameters(InToExModel);


% Create Mex-files for faster simulation
try
    SBPDmakeMEXmodel(InModel,     'InModel_MexFile');
    SBPDmakeMEXmodel(InToExModel, 'InModelToEx_MexFile');
    p.CompileFlag = 1;  % compilation successful
catch ModelCompileError
    p.CompileFlag = 0;  % compilation failed
    p.InModel = InModel;
    p.InToExModel = InToExModel;
end    
    
p.In.Vex0 = max(p.Ex.Moi, 1); % number of virus particles infecting a cell, Vex is based on MOI, but can't be lower than 1 virus particle
p.In.Ic(strcmp('Vex', p.In.StateNames)) = p.In.Vex0;                        % set initial Vex to MOI, minimum is one virus particle per cell
p.InToEx.Ic(strcmp('VpCyt', p.InToEx.StateNames)) = ...                     % set initial VpCyt in "coupling" cell to 8*MOI*Ffus
    max([8 8*p.In.Vex0*p.In.ParaValues(strcmp('Ffus', p.In.ParaNames))]);   % minimum is 8 vRNPs after infection       

% adjust Fpar(0) to MOI conditions
if ( p.Ex.Moi >= 73 )
    p.In.Ic(strcmp('FPar', p.In.StateNames))         = p.In.ChangedIcValues(strcmp('FPar',p.In.ChangedIcNames));
    p.InToEx.Ic(strcmp('FPar', p.InToEx.StateNames)) = p.In.ChangedIcValues(strcmp('FPar',p.In.ChangedIcNames));
else
    p.In.Ic(strcmp('FPar', p.In.StateNames))         = 0.26;
    p.InToEx.Ic(strcmp('FPar', p.InToEx.StateNames)) = 0.26;
end

% Simulate model with Mex-files or MatLab solver
if( p.CompileFlag )
    result.In     = SBPDsimulate('InModel_MexFile', max(d.In.Time), p.In.Ic); 
    result.InToEx = SBPDsimulate('InModelToEx_MexFile', p.Ex.tspan, p.InToEx.Ic);
else
    result.In     = SBsimulate(InModel,     'ode23s', max(d.In.Time), p.In.Ic); 
    result.InToEx = SBsimulate(InToExModel, 'ode23s', p.Ex.tspan, p.InToEx.Ic);
end

%% Prepare simulation of the extracellular model
p.Ex.ModelName = @Sim_EquationSystem;      % function handle to m-file containin the extracellular model

result.Ex.states    = {'VAttHi', 'VAttLo', 'VEn', 'T', 'Ta', 'Ia', 'V', 'Par', 'D', 'VRelTot', 'VRelEx'};
result.Ex.variables = {'I', 'Log10Tcid50', 'Log10Ha'};
p.Ex.StateNames = result.Ex.states; p.Ex.VariableNames = result.Ex.variables;

p.Idx.T      = strcmp('T',       result.Ex.states);
p.Idx.Ta     = strcmp('Ta',      result.Ex.states);
p.Idx.Ia     = strcmp('Ia',      result.Ex.states);
p.Idx.VEn    = strcmp('VEn',     result.Ex.states);
p.Idx.V      = strcmp('V',       result.Ex.states);
p.Idx.Par    = strcmp('Par',     result.Ex.states);
p.Idx.VAttHi = strcmp('VAttHi',  result.Ex.states);
p.Idx.VAttLo = strcmp('VAttLo',  result.Ex.states);
p.Idx.I      = strcmp('I',       result.Ex.variables);
p.Idx.D      = strcmp('D',       result.Ex.states);
p.Idx.VRelTot= strcmp('VRelTot', result.Ex.states);
p.Idx.VRelEx = strcmp('VRelEx',  result.Ex.states);

% initial conditions
p.Ex.y0 = zeros(1,length(result.Ex.states));
for CountStates = 1:length(d.Ex.ExpStateNames)
    p.Ex.y0(strcmpi(d.Ex.ExpStateNames(CountStates), result.Ex.states)) = ...
        d.Ex.ExpStateValues(1, CountStates);
end

% adjust initial cell + virus concentration to MOI conditions
if ( p.Ex.Moi == 73 )
    p.Ex.I0          = 1.75e5;                    % set specific initial concentration of infected cells for MOI 73 based on experimental data [Frensing, 2016]
else
    p.Ex.y0(p.Idx.T) = 5.99e5;                    % adjust total cell number for lower MOI conditions
    p.Ex.I0          = 0;                         % set initial concentration of infected cells
end

p.Ex.y0(p.Idx.V)    = p.Ex.Moi*p.Ex.y0(p.Idx.T);  % adjust starting virus concentration
p.Ex.y0(p.Idx.Par)  = 0; 

%% Simulation and Optimization
switch lower(p.Task)
    case {'simulate'}
        clc
        disp('Simulating...')

        %% Calculate the integral of kApoI(tau), exp(-kApoI(tau)) and rRel(tau)*exp(-kApoI(tau)) over tau
        p.InToEx.Rrel     = result.InToEx.reactionvalues(:,strcmp('rRel', result.InToEx.reactions));
        p.InToEx.Rrel_Par = result.InToEx.reactionvalues(:,strcmp('rRel_par', result.InToEx.reactions));

        %modify apoptosis rates and integrals
        [p.Ex.ApoI, p.Ex.IntKapoi] = getApoI(p);

        p.Ex.RrelIntKapoi     = p.InToEx.Rrel'.*p.Ex.IntKapoi;
        p.Ex.RrelIntKapoi_Par = p.InToEx.Rrel_Par'.*p.Ex.IntKapoi;

        %% Solve the system using Euler's Method
        [result.Ex.time, result.Ex.statevalues, result.Ex.variablevalues(:, p.Idx.I), ...
            result.Ex.AgeI] = Sim_Euler(p);

        result.Ex.variablevalues(:,strcmp('Log10Tcid50', result.Ex.variables)) =...
            log10(max(1, result.Ex.statevalues(:,p.Idx.V)));   
        result.Ex.variablevalues(:,strcmp('Log10Ha', result.Ex.variables)) =...
            log10(max(1, result.Ex.statevalues(:,p.Idx.Par)));          

        if ( p.Ex.Moi > 3 && p.Ex.Moi < 73 )
            fprintf('Warning: Simulation for an MOI between 3 and 73 uses parameter settings intended for\n');
            fprintf('         low/medium MOI and may not be representative for the full spectrum of MOIs.\n');
        end
        
        if ( ~strcmp('logistic',p.Apo.Increase) )
            fprintf('Warning: Simulation parameters are tuned for using a LOGISTIC FUNCTION to calculate an\n');
            fprintf('         increasing apoptosis rate. Consider optimizing the parameters when using a\n');
            fprintf('         different apoptosis increase function.\n');
        end
        
        Sim_VisualizeResults(result, d, p);

    case {'optimize'}
        clc
        if ( p.Ex.Moi ~= 73 )
            error('The parameter estimation can only be performed with an MOI of 73.');
        end
        disp('Optimizing...')
        fprintf('Start:         %s\n',datestr(p.Info.StartTime));
        if ( p.MaxTime >= 24*60*60 );  tmp = datevec(p.MaxTime./(60*60*24)) - [0 1 0 0 0 0]; 
        else                           tmp = datevec(p.MaxTime./(60*60*24));
        end
        fprintf('Projected end: %s\n',datestr(p.Info.StartTime + tmp));

        %% Optimize the multiscale model
        Opt = Opt_Optimization;
        if(isempty(Opt));  
            disp('Error!');  
            return
        end

        %% Display estimated parameter values
        if ( strcmp('fSSm', p.Solver) )
            AdjWhiteSpace = cell(length(p.problem.x_0),2);          % adjustable white space, get values in right position
            tmp = getPotency(p.problem.x_0,2);
            for i = 1 : length(p.problem.x_0)
                IdvWhiteSpace = max(tmp(:,2)) - tmp(i,2);
                if ( IdvWhiteSpace > max(tmp(:,2)) );  IdvWhiteSpace = max(tmp(:,2));  end
                for j = 1 : IdvWhiteSpace;  AdjWhiteSpace{i,1} = [AdjWhiteSpace{i,1} ' '];  end
            end
            tmp = getPotency([Opt.In.xbest, Opt.Ex.xbest],2);
            for i = 1 : length([Opt.In.xbest, Opt.Ex.xbest])
                IdvWhiteSpace = max(tmp(:,2)) - tmp(i,2);
                if ( IdvWhiteSpace > max(tmp(:,2)) );  IdvWhiteSpace = max(tmp(:,2));  end
                for j = 1 : IdvWhiteSpace;  AdjWhiteSpace{i,2} = [AdjWhiteSpace{i,2} ' '];  end
            end  
            OptNames = [p.In.PaOpt p.In.IcOpt p.Ex.PaOpt p.Ex.IcOpt];   % shorten names of parameters if they are too long
            for i = 1 : length(tmp)
                if ( length(OptNames{i}) > 7 );  OptNames{i} = OptNames{i}(1:7);  end
            end

            fprintf('\n+++++++++++ Parameter estimation +++++++++++\n');
            fprintf('----initial----        ----estimated----\n');
            aws = 1;
            for i = 1 : size(p.In.PaOpt,2);  fprintf('%s%.6f\t%s\t%s%.9f\n',AdjWhiteSpace{aws,1},p.problem.x_0(aws),OptNames{aws},AdjWhiteSpace{aws,2},Opt.In.xbest(i));  aws = aws + 1; end
            for i = 1 : size(p.In.IcOpt,2);  fprintf('%s%.6f\t%s\t%s%.9f\n',AdjWhiteSpace{aws,1},p.problem.x_0(aws),OptNames{aws},AdjWhiteSpace{aws,2},Opt.In.xbest(size(p.In.PaOpt,2)+i));  aws = aws + 1; end
            fprintf('++++++++++++++++++++++++++++++++++++++++++++\n');
            for i = 1 : size(p.Ex.PaOpt,2);  fprintf('%s%.6f\t%s\t%s%.9f\n',AdjWhiteSpace{aws,1},p.problem.x_0(aws),OptNames{aws},AdjWhiteSpace{aws,2},Opt.Ex.xbest(i));  aws = aws + 1;  end        
            for i = 1 : size(p.Ex.IcOpt,2);  fprintf('%s%.6f\t%s\t%s%.9f\n',AdjWhiteSpace{aws,1},p.problem.x_0(aws),OptNames{aws},AdjWhiteSpace{aws,2},Opt.Ex.xbest(size(p.Ex.PaOpt,2)+i));  aws = aws + 1;  end  
            fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');
        end

        %% Simulate the optimized intracellular model
        p.In.Ic(p.In.PosIcOpt) = Opt.In.xbest(p.In.nPaOpt+1:end);
        if ( p.InToEx.PosIcOpt(1) == 1 )
            p.InToEx.Ic(1) = max([8 8*Opt.In.xbest(p.In.nPaOpt+1)*p.In.ParaValues(strcmp('Ffus', p.In.ParaNames))]);
            if ( length(p.InToEx.PosIcOpt) > 2 )
                p.InToEx.Ic(p.InToEx.PosIcOpt(2:end)) = Opt.In.xbest(p.In.nPaOpt+2:p.In.nOpt);
            end
        else        
            p.InToEx.Ic(p.InToEx.PosIcOpt) = Opt.In.xbest(p.In.nPaOpt+1:p.In.nOpt);         
        end

        clear mex
        h.In.PaBest = Opt.In.xbest(1:p.In.nPaOpt);

        % Simulate model with Mex-files or MatLab solver
        if( p.CompileFlag )
            result.In     = SBPDsimulate('InModel_MexFile', max(d.In.Time)*1.05, p.In.Ic, p.In.PaOpt, h.In.PaBest); 
            result.InToEx = SBPDsimulate('InModelToEx_MexFile', p.Ex.tspan, p.InToEx.Ic, p.In.PaOpt(~p.In.IdxKfus), h.In.PaBest(~p.In.IdxKfus));
        else
            InModel_tmp = SBparameters(p.In.PaOpt, h.In.PaBest);
            InToExModel_tmp = SBparameters(p.In.PaOpt(~p.In.IdxKfus), h.In.PaBest(~p.In.IdxKfus));
            result.In     = SBsimulate(InModel,     'ode23s', max(d.In.Time)*1.05, p.In.Ic); 
            result.InToEx = SBsimulate(InToExModel, 'ode23s', p.Ex.tspan, p.InToEx.Ic);
        end        
        
        p.InToEx.Rrel     = result.InToEx.reactionvalues(:,strcmp('rRel', result.InToEx.reactions));
        p.InToEx.Rrel_Par = result.InToEx.reactionvalues(:,strcmp('rRel_par', result.InToEx.reactions));

        %% Simulate the optimized extracellular model
        p.Ex.y0(p.Ex.PosIcOpt) = Opt.Ex.xbest(p.Ex.nPaOpt+1:end);
        for CountPara = 1:p.Ex.nPaOpt
            p.Ex.(char(p.Ex.PaOpt(CountPara))) = Opt.Ex.xbest(CountPara);
        end
        if(sum(strcmp('kFus', p.In.PaOpt)));   p.Ex.kFus  = Opt.In.xbest(p.In.IdxKfus);   end
        if(sum(strcmp('kApoT', p.In.PaOpt)));  p.Ex.kApoT = Opt.In.xbest(p.In.IdxKapoT);  end
        if(sum(strcmp('kApoI', p.In.PaOpt)));  p.Ex.kApoI = Opt.In.xbest(p.In.IdxKapoI);  end
        if(sum(strcmp('ASS', p.In.PaOpt)));    p.Ex.ASS   = Opt.In.xbest(p.In.IdxASS);    end
        if(sum(strcmp('AST', p.In.PaOpt)));    p.Ex.AST   = Opt.In.xbest(p.In.IdxAST);    end

        %modify apoptosis rates and integrals
        [p.Ex.ApoI, p.Ex.IntKapoi] = getApoI(p);

        p.Ex.RrelIntKapoi     = p.InToEx.Rrel'.*p.Ex.IntKapoi;
        p.Ex.RrelIntKapoi_Par = p.InToEx.Rrel_Par'.*p.Ex.IntKapoi;

        %simulate the system with new parameters and initial conditions
        [result.Ex.time, result.Ex.statevalues, result.Ex.variablevalues(:, p.Idx.I)]...
            = Sim_Euler(p);

        result.Ex.variablevalues(:,strcmp('Log10Tcid50', result.Ex.variables)) =...
            log10(max(1, result.Ex.statevalues(:,p.Idx.V)));
        result.Ex.variablevalues(:,strcmp('Log10Ha', result.Ex.variables)) =...
            log10(max(1, result.Ex.statevalues(:,p.Idx.Par)));         

        result.Opt = Opt;

        Sim_VisualizeResults(result, d, p);

        %% Display estimation information
        p.Info.EndTime = clock;
        fprintf('\n++++++++++ Additional information ++++++++++\n');
        fprintf('Start:      %s\n',datestr(p.Info.StartTime));
        fprintf('End:        %s\n',datestr(p.Info.EndTime));
        if ( datenum(p.Info.EndTime - p.Info.StartTime) < 0.0417 )
            fprintf('Duration:   %i min\n',ceil(1440*datenum(p.Info.EndTime - p.Info.StartTime)));
        else
            tmp = 24*datenum(p.Info.EndTime - p.Info.StartTime);
            fprintf('Duration:   %i h + %i min\n',floor(tmp),ceil(60*(tmp-floor(tmp))));
        end
        try
            fprintf('FunEvals:   %i\n',result.Opt.fSSm.numeval);
            fprintf('ObjFunVal:  %i\n',result.Opt.ObjFunVal/1e6);
        catch ME
        end
        fprintf('++++++++++++++++++++++++++++++++++++++++++++\n\n');

    otherwise
        disp('Warnin(Main): No such task.')
end
   
%% Housekeeping commands II
clear Count* h
path(p.SearchPath)

































