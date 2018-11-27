function Sim_VisualizeResults(result, d, p)
%SIM_VISUALIZERESULTS plots the simulation results and adds the
%   experimental data if available for the specified MOI
%
%   last revised: 2018/11/21

%% Retrieve states
SimTimeIn   = result.In.time;
SimStatesIn = [result.In.states result.In.variables];
SimValuesIn = [result.In.statevalues result.In.variablevalues];

SimTimeEx   = result.Ex.time;
SimStatesEx = [result.Ex.states result.Ex.variables];
SimValuesEx = [result.Ex.statevalues result.Ex.variablevalues];

SimIndex.Mrna    = strcmp('Rm5',         SimStatesIn);
SimIndex.Crna    = strcmp('RcSegment',   SimStatesIn);
SimIndex.Vrna    = strcmp('RvSegment',   SimStatesIn);
SimIndex.Vrel    = strcmp('Vrel',        SimStatesIn);
SimIndex.Prel    = strcmp('Prel',        SimStatesIn);
SimIndex.FPar    = strcmp('FPar',        SimStatesIn);
SimIndex.T       = strcmp('T',           SimStatesEx);
SimIndex.I       = strcmp('I',           SimStatesEx);
SimIndex.Ta      = strcmp('Ta',          SimStatesEx);
SimIndex.Ia      = strcmp('Ia',          SimStatesEx);
SimIndex.V       = strcmp('V',           SimStatesEx);
SimIndex.Par     = strcmp('Par',         SimStatesEx);
SimIndex.VRelTot = strcmp('VRelTot',     SimStatesEx);
SimIndex.VEn     = strcmp('VEn',         SimStatesEx);
SimIndex.Tcid50  = strcmp('Log10Tcid50', SimStatesEx);
SimIndex.log10Ha = strcmp('Log10Ha',     SimStatesEx);
SimIndex.VAttHi  = strcmp('VAttHi',      SimStatesEx);
SimIndex.VAttLo  = strcmp('VAttLo',      SimStatesEx);
SimIndex.D      = strcmp('D',           SimStatesEx);

if ( p.Ex.Moi == 73 )
    ExpIndex.Mrna   = strcmp('Rm5',          d.In.ExpStateNames);
    ExpIndex.Crna   = strcmp('RcSegment',    d.In.ExpStateNames);
    ExpIndex.Vrna   = strcmp('RvSegment',    d.In.ExpStateNames);
    ExpIndex.Vrel   = strcmp('Vrel',         d.In.ExpStateNames);
    ExpIndex.Prel   = strcmp('Prel',         d.In.ExpStateNames);
    ExpIndex.FPar   = strcmp('FPar',         d.In.ExpStateNames);
    ExpIndex.VrnaSd = strcmp('StdRvSegment', d.In.ExpStateNames);
    ExpIndex.CrnaSd = strcmp('StdRcSegment', d.In.ExpStateNames);
    ExpIndex.MrnaSd = strcmp('StdRm5',       d.In.ExpStateNames);
    ExpIndex.VrelSd = strcmp('StdVrel',      d.In.ExpStateNames);
    ExpIndex.PrelSd = strcmp('StdPrel',      d.In.ExpStateNames);
    ExpIndex.FParSd = strcmp('StdFPar',      d.In.ExpStateNames);
end    
ExpIndex.T      = strcmp('T',            d.Ex.ExpStateNames);
ExpIndex.I      = strcmp('I',            d.Ex.ExpStateNames);
ExpIndex.Ia     = strcmp('Ia',           d.Ex.ExpStateNames);
ExpIndex.Ta     = strcmp('Ta',           d.Ex.ExpStateNames);
ExpIndex.V      = strcmp('V',            d.Ex.ExpStateNames);
ExpIndex.Par    = strcmp('Par',          d.Ex.ExpStateNames);
ExpIndex.TSd    = strcmp('StdT',         d.Ex.ExpStateNames);
ExpIndex.ISd    = strcmp('StdI',         d.Ex.ExpStateNames);
ExpIndex.IaSd   = strcmp('StdIa',        d.Ex.ExpStateNames);
ExpIndex.TaSd   = strcmp('StdTa',        d.Ex.ExpStateNames);
ExpIndex.VSd    = strcmp('StdV',         d.Ex.ExpStateNames);
ExpIndex.ParSd  = strcmp('StdPar',       d.Ex.ExpStateNames);



%% realize offset in RNA values
% as free viral RNAs in the seed virus might adhere to cells and thus cause
% a constant offset in measurements we implement such an offset
if ( p.Ex.Moi == 73 )
    RnaNames = {'Rm5', 'RcSegment', 'RvSegment'};
    p.In.OffsetIndex = zeros(length(RnaNames), 2);
    for CountRnas = 1:length(RnaNames)
        p.In.OffsetIndex(CountRnas, 1) = find(strcmp(RnaNames(CountRnas), SimStatesIn));
        p.In.OffsetIndex(CountRnas, 2) = find(strcmp(RnaNames(CountRnas), d.In.ExpStateNames));
    end

    SimValuesIn(:,p.In.OffsetIndex(:,1)) = SimValuesIn(:,p.In.OffsetIndex(:,1)) +...
        repmat(d.In.ExpStateValues(1,p.In.OffsetIndex(:,2)), length(SimTimeIn), 1);
    
    disp('Warning: Offset of first data point for all RNA values implemented.');
end

%% Visualize - vRNA, cRNA, mRNA and released virions
h.Fig = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [1 1 22 10],...
               'Name', 'Intracellular Level');
                     
%% mRNA level
h.Ax1 = subplot(231);
set(h.Ax1, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.08 .6 .22 .32]);

plot(h.Ax1, SimTimeIn, SimValuesIn(:,SimIndex.Mrna), 'b', 'LineWidth', 1.5); hold on
if ( p.Ex.Moi == 73 )
    errorbar(d.In.Time, d.In.ExpStateValues(:,ExpIndex.Mrna), d.In.ExpStateValues(:,ExpIndex.MrnaSd), 'bo', 'LineWidth', 1, 'MarkerSize', 6); hold off
end

set(h.Ax1, 'Xlim', [0 max(SimTimeIn)*1.05], 'XTick', 0:ceil(max(SimTimeIn)/6):max(SimTimeIn),...
           'YLim', [1e-1 1e5], 'YTick', 10.^(-2:5),'YScale', 'log');
       
xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('molecules per cell', 'FontName', 'Arial', 'FontSize', 10);

text(9, 3e2, 'mRNA', 'Color', 'b',...
    'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center');

%% cRNA level
h.Ax2 = subplot(232);
set(h.Ax2, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.4 .6 .22 .32]);

plot(h.Ax2, SimTimeIn, SimValuesIn(:,SimIndex.Crna), 'r', 'LineWidth', 1.5); hold on;
if ( p.Ex.Moi == 73 )
    errorbar(d.In.Time, d.In.ExpStateValues(:,ExpIndex.Crna), d.In.ExpStateValues(:,ExpIndex.CrnaSd), 'ro', 'LineWidth', 1, 'MarkerSize', 6); hold off
end

set(h.Ax2, 'Xlim', [0 max(SimTimeIn)*1.05], 'XTick', 0:ceil(max(SimTimeIn)/6):max(SimTimeIn),...
           'YLim', [1e-1 1e3], 'YTick', 10.^(-2:5),'YScale', 'log');     
       
xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('molecules per cell', 'FontName', 'Arial', 'FontSize', 10);

text(10, 5e1, 'cRNA', 'Color', 'r',...
    'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center');

%% vRNA level
h.Ax3 = subplot(233);
set(h.Ax3, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.72 .6 .22 .32]);
       
plot(h.Ax3, SimTimeIn, SimValuesIn(:,SimIndex.Vrna), 'k', 'LineWidth', 1.5); hold on;
if ( p.Ex.Moi == 73 )
    errorbar(d.In.Time, d.In.ExpStateValues(:,ExpIndex.Vrna), d.In.ExpStateValues(:,ExpIndex.VrnaSd), 'ko', 'LineWidth', 1, 'MarkerSize', 6); hold off
end


set(h.Ax3, 'Xlim', [0 max(SimTimeIn)*1.05], 'XTick', 0:ceil(max(SimTimeIn)/6):max(SimTimeIn),...
           'YLim', [1e1 1e5], 'YTick', 10.^(-2:5),'YScale', 'log');     
       
xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('molecules per cell', 'FontName', 'Arial', 'FontSize', 10);

text(10, 2e2, 'vRNA', 'Color', 'k',...
    'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center');

%% released virions

%%% infective virus particles
h.Ax4 = subplot(234);
set(h.Ax4, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.22 .12 .22 .32]);

plot(h.Ax4, SimTimeIn, SimValuesIn(:,SimIndex.Vrel), 'Color', [0 .5 0], 'LineWidth', 1.5); hold on;
if ( p.Ex.Moi == 73 )
    errorbar(d.In.Time, d.In.ExpStateValues(:,ExpIndex.Vrel), d.In.ExpStateValues(:,ExpIndex.VrelSd), 'o', 'Color', [0 .5 0], 'LineWidth', 1, 'MarkerSize', 6); hold off
    ylim([0 max([max(SimValuesIn(:,SimIndex.Vrel)), ...
                   max(d.In.ExpStateValues(:,ExpIndex.Vrel)+d.In.ExpStateValues(:,ExpIndex.VrelSd))])*1.05])
else
    ylim([0 max(SimValuesIn(:,SimIndex.Vrel))*1.05])
end

set(h.Ax4, 'Xlim', [0 max(SimTimeIn)*1.05], 'XTick', 0:ceil(max(SimTimeIn)/6):max(SimTimeIn),...
           'YTick', 0:100*ceil(max(SimValuesIn(:,SimIndex.Vrel))/600):max(SimValuesIn(:,SimIndex.Vrel))*10);

       
xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('virions per cell', 'FontName', 'Arial', 'FontSize', 10);

text(7, 0.7*max(SimValuesIn(:,SimIndex.Vrel)), sprintf('infectious\nparticles'), 'Color', [0 .5 0],...
    'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center');

%%% non-infective virus particles
h.Ax5 = subplot(236);
set(h.Ax5, 'FontSize', 9,...
            'FontName', 'Arial',...
            'LineWidth', 1,...
            'Position', [.56 .12 .22 .32]);

plot(h.Ax5, SimTimeIn, SimValuesIn(:,SimIndex.Prel), 'Color', [0 .5 0], 'LineWidth', 1.5); hold on;
if ( p.Ex.Moi == 73 )
    errorbar(d.In.Time, d.In.ExpStateValues(:,ExpIndex.Prel), d.In.ExpStateValues(:,ExpIndex.PrelSd), 'o', 'Color', [0 .5 0], 'LineWidth', 1, 'MarkerSize', 6); hold off
    ylim([0 max([max(SimValuesIn(:,SimIndex.Prel)),...
                           max(d.In.ExpStateValues(:,ExpIndex.Prel)+d.In.ExpStateValues(:,ExpIndex.PrelSd))])*1.05])
else
    ylim([0 max(SimValuesIn(:,SimIndex.Prel))*1.05]);
end

set(h.Ax5, 'Xlim', [0 max(SimTimeIn)*1.05], 'XTick', 0:ceil(max(SimTimeIn)/6):max(SimTimeIn),...
           'YTick', 0:5e3:5e4);

xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('virions per cell', 'FontName', 'Arial', 'FontSize', 10);

text(7, 0.7*max(SimValuesIn(:,SimIndex.Prel)), sprintf('total\nparticles'), 'Color', [0 .5 0],...
    'FontName', 'Arial', 'FontSize', 12, 'HorizontalAlignment', 'center');     

%% Visualization - infection dynamics and virus titers
h.Fig2 = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 5],...
               'units', 'centimeters',...
               'position', [1 14 22 8],...
               'Name', 'Extracellular Level');         

h.Ax5 = subplot(221);
set(h.Ax5, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.08 .15 .48 .7]);
       
h.L1 = plot(h.Ax5, SimTimeEx, SimValuesEx(:,SimIndex.T),  'b',...
                   SimTimeEx, SimValuesEx(:,SimIndex.I),  'g',...
                   SimTimeEx, SimValuesEx(:,SimIndex.Ia), 'k-',...
                   'LineWidth', 1.5); set(h.L1(2), 'Color', [0 .5 0]); hold on;                   
if ( p.Ex.Moi == 73 ) 
    h.L2 = plot(h.Ax5, d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.T), 'bo',...
                       d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.I), 'go',...
                       d.Ex.Time, d.Ex.ExpStateValues(:, ExpIndex.Ia), 'ko',...
                       'LineWidth', 1.5); set(h.L2(2), 'Color', [0 .5 0]); 
end

set(h.Ax5, 'XLim', [0 max(SimTimeEx)], 'XTick', 0:12:max(SimTimeEx),...           
           'YLim', [0 max([max(SimValuesEx(:,SimIndex.T|SimIndex.I|SimIndex.Ta|SimIndex.Ia)), max(d.Ex.ExpStateValues(:, ExpIndex.T))])*1.05]);

xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('cells per ml', 'FontName', 'Arial', 'FontSize', 10); 
title('    infection dynamics', 'FontName', 'Arial', 'FontSize', 12);
h.Leg = legend('target cells', 'infected cells', 'apo. infected cells', 'Location','EastOutside');
set(h.Leg, 'FontName', 'Arial', 'FontSize', 10)

%%%
h.Ax6 = subplot(222);
set(h.Ax6, 'FontSize', 9,...
           'FontName', 'Arial',...
           'LineWidth', 1,...
           'Position', [.65 .15 .27 .7]);      
       
% TCID50       
line_TCID50 = plot(h.Ax6, SimTimeEx, log10(max(1, result.Ex.statevalues(:,p.Idx.VRelTot))),...
                          'b','LineWidth', 1.5); hold on;

if ( sum(p.Ex.Moi == [73 1e-4]) )
    TCID50_lb = zeros(1,length(d.Ex.ExpStateValues(:, ExpIndex.V)));
    TCID50_ub = zeros(1,length(d.Ex.ExpStateValues(:, ExpIndex.V)));


    % calculate lower + upper bounds for TCID50
    TCID50_lb = log10(d.Ex.ExpStateValues(:, ExpIndex.V)) - log10(d.Ex.ExpStateValues(:,ExpIndex.V) - d.Ex.ExpStateValues(:,ExpIndex.VSd));
    for i = 1 : length(TCID50_lb)
        if ( ~isreal(TCID50_lb(i)) || isinf(TCID50_lb(i)) );  TCID50_lb(i) = log10(d.Ex.ExpStateValues(i, ExpIndex.V));  end
    end

    TCID50_ub = log10(d.Ex.ExpStateValues(:, ExpIndex.V) + d.Ex.ExpStateValues(:,ExpIndex.VSd)) - log10(d.Ex.ExpStateValues(:, ExpIndex.V));



    errorbar(d.Ex.Time(2:end), log10(d.Ex.ExpStateValues(2:end, ExpIndex.V)),...
                        TCID50_lb(2:end), ...
                        TCID50_ub(2:end), 'bo', 'LineWidth', 1, 'MarkerSize', 6);   
elseif ( p.Ex.Moi == 3 )
    h.p = plot(d.Ex.Time, log10(d.Ex.ExpStateValues(:,1:3)),'bo','MarkerSize',6);
    set(h.p(2),'Marker','x');
    set(h.p(3),'Marker','d');
else
end
                                
%HA 
line_HA = plot(h.Ax6, SimTimeEx, SimValuesEx(:,SimIndex.log10Ha), 'r','LineWidth', 1.5);

if ( sum(p.Ex.Moi == [73 1e-4]) )
    % calculate lower + upper bounds for TCID50
    HA_lb = log10(d.Ex.ExpStateValues(:, ExpIndex.Par)) - log10(d.Ex.ExpStateValues(:,ExpIndex.Par) - d.Ex.ExpStateValues(:,ExpIndex.ParSd));
    for i = 1 : length(HA_lb)
        if ( ~isreal(HA_lb(i)) || isinf(HA_lb(i)) );  HA_lb(i) = log10(d.Ex.ExpStateValues(i, ExpIndex.Par));  end
    end

    HA_ub = log10(d.Ex.ExpStateValues(:, ExpIndex.Par) + d.Ex.ExpStateValues(:,ExpIndex.ParSd)) - log10(d.Ex.ExpStateValues(:, ExpIndex.Par));
    errorbar(d.Ex.Time, log10(d.Ex.ExpStateValues(:, ExpIndex.Par)),...
                        HA_lb, ...
                        HA_ub, 'ro', 'LineWidth', 1, 'MarkerSize', 6); 
elseif ( p.Ex.Moi == 3 )
    h.p = plot(d.Ex.Time, log10(d.Ex.ExpStateValues(:,4:6)),'ro','MarkerSize',6);
    set(h.p(2),'Marker','x');
    set(h.p(3),'Marker','d');
else
end

set(h.Ax6, 'XLim', [0 max(SimTimeEx)], 'XTick', 0:12:max(SimTimeEx),... % );
           'YLim', [3.5 max(SimValuesEx(:,SimIndex.log10Ha))*1.05]);
       
if ( p.Ex.Moi == 1e-4 )
    ylim([1 max(SimValuesEx(:,SimIndex.log10Ha))*1.05])
end
       
xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel('log_{10} particles per ml', 'FontName', 'Arial', 'FontSize', 10); 
title('virus titers', 'FontName', 'Arial', 'FontSize', 12);
h.Leg2 = legend([line_TCID50 line_HA],'infectious', 'total',  'Location', 'SouthEast'); % only show legend for simulation
set(h.Leg2, 'FontName', 'Arial', 'FontSize', 10)




