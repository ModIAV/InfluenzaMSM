function para = Sim_ParameterDeclaration
%SIM_PARAMETERDECLARATION defines and saves parameters for an age-segregated 
%   model of influenza A virus infection
%
%   last revised: 2018/11/21

%% Attachment and internalization parameters
para.Finf = 1;              %cells/virion  number of cells that are infected by one infectious virion in the cytoplasm (equal or lower than 1)

para.BtotHi  = 150;         %sites/cell  total amount of high-affinity binding sites on the cell surface (according to Nunes-Correia et al. 1999)
para.BtotLo  = 1000;        %sites/cell  total amount of low-affinity binding sites on the cell surface (according to Nunes-Correia et al. 1999)
para.kAttcHi = 3.32e-8;     %mL/(site*h)  attachment rate of virus particles to high-affinity binding sites (Heldt et al. 2012)
para.kAttcLo = 1.85e-10;    %mL/(site*h)  attachment rate of virus particles to low-affinity binding sites (Heldt et al. 2012)
para.kEqcHi  = 4.48e-9;     %mL/site  equilibrium constant for the attachment of virus particles to the high-affinitiy binding sites (Nunes-Correia et al. 1999)
para.kEqcLo  = 3.32e-11;    %mL/site  equilibrium constant for the attachment of virus particles to the low-affinitiy binding sites (Nunes-Correia et al. 1999)
para.kEn     = 4.8;         %1/h  endocytosis rate or virions bound to the high-affinity and low-affinity binding sites
para.kFus    = 9.56e-3;     %1/h  fusion rate of virions in late endosomes

para.kDisHi = para.kAttcHi/para.kEqcHi; %1/h  dissociaten rate of virus bound to high-affinity binding sites
para.kDisLo = para.kAttcLo/para.kEqcLo; %1/h  dissociaten rate of virus bound to low-affinity binding sites

%% Growth and infection dynamics parameters
para.Tmax  = 10e5;         %cells/ml  maximum cell concentration supported by the surface of the culture vessel

para.MuMax = 0.03;         %1/h  maximum growth rate of target cells (Schulze-Horsel et al. 2009, RKI)
para.kApoT = 7.35e-3;      %1/h  apoptosis rate of target cells
para.kApoI = 3.28e-2;      %1/h  increase in apoptosis rate caused by viral infection
para.kLys  = 6.39e-2;      %1/h  death/lysis rate of apoptotic cells 
para.kDegV = 0.1;          %1/h  degradation rate of infectious virions (Beauchemin 2008)

%% New parameters
para.v_Apo  = 1;            % variance of infected cell lifespan
para.tau_Apo  = 20;         % average infected cell lifespan

para.ParaNames = fieldnames(para);


