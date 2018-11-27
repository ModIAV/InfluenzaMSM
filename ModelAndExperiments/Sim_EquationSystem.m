function dydt = Sim_EquationSystem(~, y, IntRrelI, IntRrelI_Par, IntI, IntApoI, p)
%SIM_EQUATIONSYSTEM contains the model equations of the age-segregated 
%   model for influenza virus infection at the extracellular level
%
%   last revised: 2018/11/21

%% retrieve state variables
T       = y(p.Idx.T);
Ta      = y(p.Idx.Ta);
Ia      = y(p.Idx.Ia);
V       = y(p.Idx.V);
Par     = y(p.Idx.Par);
VEn     = y(p.Idx.VEn);
VAttHi  = y(p.Idx.VAttHi);
VAttLo  = y(p.Idx.VAttLo);
D       = y(p.Idx.D);
VRelTot = y(p.Idx.VRelTot);
VRelEx = y(p.Idx.VRelEx);

%% algebraic equations
mu  = max(0, p.Ex.MuMax/p.Ex.Tmax * (p.Ex.Tmax - (T + IntI)));
BHi = max(0, p.Ex.BtotHi*(T+Ta) - VAttHi);    %sites/ml  concentration of free high-affinity binding sites
BLo = max(0, p.Ex.BtotLo*(T+Ta) - VAttLo);    %sites/ml  concentration of free low-affinity binding sites

if(T+Ta>0)
    rInf = p.Ex.Finf*p.Ex.kFus*VEn/(T+Ta);
    rLys = p.Ex.kLys*Ta           /(T+Ta);
else
    rInf = 0; rLys = 0;
end

%% differential equations
dydt = zeros(size(y));

%attachment to the cell surface
dydt(p.Idx.VAttHi)  = p.Ex.kAttcHi*BHi*V - (p.Ex.kDisHi + p.Ex.kEn)*VAttHi - (rInf + rLys)*VAttHi;
dydt(p.Idx.VAttLo)  = p.Ex.kAttcLo*BLo*V - (p.Ex.kDisLo + p.Ex.kEn)*VAttLo - (rInf + rLys)*VAttLo;

%endocytosis and fusion
dydt(p.Idx.VEn)     = p.Ex.kEn*(VAttHi+VAttLo) - p.Ex.kFus*VEn             - (rInf + rLys)*VEn;

%cell populations
dydt(p.Idx.T)       = (mu          - rInf    - p.Ex.kApoT)*T;
dydt(p.Idx.Ta)      = p.Ex.kApoT*T - rInf*Ta - p.Ex.kLys  *Ta;
dydt(p.Idx.Ia)      = IntApoI      + rInf*Ta - p.Ex.kLys  *Ia;

%virus concentrations
dydt(p.Idx.V)       = IntRrelI - p.Ex.kDegV*V + p.Ex.kDisHi*VAttHi + p.Ex.kDisLo*VAttLo - (p.Ex.kAttcHi*BHi + p.Ex.kAttcLo*BLo)*V;
dydt(p.Idx.Par)     = IntRrelI_Par + p.Ex.kDegV*V;

% new 
dydt(p.Idx.D)       = p.Ex.kLys  *Ta + p.Ex.kLys  *Ia;
dydt(p.Idx.VRelTot) = IntRrelI - p.Ex.kDegV*VRelTot;

dydt(p.Idx.VRelEx)  = IntRrelI;
