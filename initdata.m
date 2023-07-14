%===variable descriptions===%

% parameters and counters
% alpha       % capital-share parameter in production function
% delta       % depreciation rate
% sigL        % elasticity of hours-disutility
% fag         % first economically active age group
% rho         % discount rate households
% sigma       % elasticity of intertemporal substitution
% tt          % time index

% age-dependent variables
% parlv0      % multiplicative shift parameter hours-disutillity
% parlv1      % additive shift parameter hours-disutillity

% time-dependent variables
% A          % aggregate assets
% Cons       % aggregate consumption
% CG         % aggregate public consumption
% DG         % public debt
% Div        % dividends
% Exp        % government expenditure
% Inv        % investment
% K          % capital stock
% LD         % aggregate labor demand
% LS         % aggregate labor supply
% N          % aggregate population
% NB         % number of newborns
% Nc         % number of children
% Nr         % number of retirees
% Nw         % number of workers
% P          % pension expenditure
% PB         % primary balance
% Rev        % government revenue
% TaxF       % payroll tax revenue
% TFP        % TFP stock
% V          % value of representative firm
% Y          % output (GDP)
% ab         % aggregate accidental bequests
% eda        % excess demand asset market
% edab       % (excess demand) accidental bequests
% edg        % (excess demand) government budget
% ediv       % (excess demand) intervivo transfers
% edl        % excess demand labor market
% edw        % (excess demand) Walras' Law
% edy        % excess demand goods market
% iv         % aggregate received intervivo transfers
% pc         % price of consumption
% qTob       % Tobin's q
% r          % real interest rate
% rag        % retirement age group
% tauC       % consumption tax rate
% tauF       % payroll tax rate
% tauW       % average wage tax rate
% taul       % lump-sum tax rate
% tauprof    % corporate tax rate
% uck        % user cost of capital
% w          % wage rate

% age-dependent and time-dependent variables (with v and z suffix)
% A          % assets
% Cons       % consumption
% HH_nonconv % counter for non-converging households
% N          % population mass
% Sav        % end-of-period savings
% ab         % accidental bequest
% dis_tot    % disutility of labor
% cG         % public consumption
% ell        % hours worked
% gam        % conditional survival probability
% iv         % intervivo transfers
% lambda     % shadow price of assets
% notret     % not retired indicator
% p          % pension income
% pc         % price of consumption
% r          % real interest rate
% tauC       % consumption tax rate
% tauW       % wage tax rate
% taul       % lump-sum tax rate
% theta      % productivity (age-specific)
% w          % wage rate
% y          % per-period labor and pension income

%=====initialization=====%
 
global ncoh; ncoh = tend + (nag-1);

% parameters and counters
global alpha;                 alpha                    = 0.0;
global delta;                 delta                    = 0.0;
global sigL;                  sigL                     = 0.0;
global fag;                   fag                      = 0.0;
global rho;                   rho                      = 0.0;
global sigma;                 sigma                    = 0.0;
global tt;                    tt                       = 0.0;

% age-dependent variables
global parlv0;                parlv0                   = zeros(nag,1);
global parlv1;                parlv1                   = zeros(nag,1);

% time-dependent variables
global A;                     A                        = zeros(1,tend);
global A0;                    A0                       = 0.0;
global Cons;                  Cons                     = zeros(1,tend);
global Cons0;                 Cons0                    = 0.0;
global CG;                    CG                       = zeros(1,tend);
global CG0;                   CG0                      = 0.0;
global DG;                    DG                       = zeros(1,tend);
global DG0;                   DG0                      = 0.0;
global Div;                   Div                      = zeros(1,tend);
global Div0;                  Div0                     = 0.0;
global Exp;                   Exp                      = zeros(1,tend);
global Exp0;                  Exp0                     = 0.0;
global Inv;                   Inv                      = zeros(1,tend);
global Inv0;                  Inv0                     = 0.0;
global K;                     K                        = zeros(1,tend);
global K0;                    K0                       = 0.0;
global LD;                    LD                       = zeros(1,tend);
global LD0;                   LD0                      = 0.0;
global LS;                    LS                       = zeros(1,tend);
global LS0;                   LS0                      = 0.0;
global N;                     N                        = zeros(1,tend);
global N0;                    N0                       = 0.0;
global NB;                    NB                       = zeros(1,tend);
global NB0;                   NB0                      = 0.0;
global Nc;                    Nc                       = zeros(1,tend);
global Nc0;                   Nc0                      = 0.0;
global Nr;                    Nr                       = zeros(1,tend);
global Nr0;                   Nr0                      = 0.0;
global Nw;                    Nw                       = zeros(1,tend);
global Nw0;                   Nw0                      = 0.0;
global P;                     P                        = zeros(1,tend);
global P0;                    P0                       = 0.0;
global PB;                    PB                       = zeros(1,tend);
global PB0;                   PB0                      = 0.0;
global Rev;                   Rev                      = zeros(1,tend);
global Rev0;                  Rev0                     = 0.0;
global TaxF;                  TaxF                     = zeros(1,tend);
global TaxF0;                 TaxF0                    = 0.0;
global TFP;                   TFP                      = zeros(1,tend);
global TFP0;                  TFP0                     = 0.0;
global V;                     V                        = zeros(1,tend);
global V0;                    V0                       = 0.0;
global Y;                     Y                        = zeros(1,tend);
global Y0;                    Y0                       = 0.0;
global ab;                    ab                       = zeros(1,tend);
global ab0;                   ab0                      = 0.0;
global eda;                   eda                      = zeros(1,tend);
global eda0;                  eda0                     = 0.0;
global edab;                  edab                     = zeros(1,tend);
global edab0;                 edab0                    = 0.0;
global edg;                   edg                      = zeros(1,tend);
global edg0;                  edg0                     = 0.0;
global ediv;                  ediv                     = zeros(1,tend);
global ediv0;                 ediv0                    = 0.0;
global edl;                   edl                      = zeros(1,tend);
global edl0;                  edl0                     = 0.0;
global edw;                   edw                      = zeros(1,tend);
global edw0;                  edw0                     = 0.0;
global edy;                   edy                      = zeros(1,tend);
global edy0;                  edy0                     = 0.0;
global iv;                    iv                       = zeros(1,tend);
global iv0;                   iv0                      = 0.0;
global pc;                    pc                       = zeros(1,tend);
global pc0;                   pc0                      = 0.0;
global qTob;                  qTob                     = zeros(1,tend);
global qTob0;                 qTob0                    = 0.0;
global r;                     r                        = zeros(1,tend);
global r0;                    r0                       = 0.0;
global rag;                   rag                      = zeros(1,tend);
global rag0;                  rag0                     = 0.0;
global tauC;                  tauC                     = zeros(1,tend);
global tauC0;                 tauC0                    = 0.0;
global tauF;                  tauF                     = zeros(1,tend);
global tauF0;                 tauF0                    = 0.0;
global tauW;                  tauW                     = zeros(1,tend);
global tauW0;                 tauW0                    = 0.0;
global taul;                  taul                     = zeros(1,tend);
global taul0;                 taul0                    = 0.0;
global tauprof;               tauprof                  = zeros(1,tend);
global tauprof0;              tauprof0                 = 0.0;
global uck;                   uck                      = zeros(1,tend);
global uck0;                  uck0                     = 0.0;
global w;                     w                        = zeros(1,tend);
global w0;                    w0                       = 0.0;

% age-dependent and time-dependent variables
global Av;                    Av                       = zeros(nag,tend);
global Az;                    Az                       = zeros(nag,ncoh);
global Av0;                   Av0                      = zeros(nag,1);
global Consv;                 Consv                    = zeros(nag,tend);
global Consz;                 Consz                    = zeros(nag,ncoh);
global Consv0;                Consv0                   = zeros(nag,1);
global HH_nonconvv;           HH_nonconvv              = zeros(nag,tend);
global HH_nonconvz;           HH_nonconvz              = zeros(nag,ncoh);
global HH_nonconvv0;          HH_nonconvv0             = zeros(nag,1);
global Nv;                    Nv                       = zeros(nag,tend);
global Nz;                    Nz                       = zeros(nag,ncoh);
global Nv0;                   Nv0                      = zeros(nag,1);
global Savv;                  Savv                     = zeros(nag,tend);
global Savz;                  Savz                     = zeros(nag,ncoh);
global Savv0;                 Savv0                    = zeros(nag,1);
global abv;                   abv                      = zeros(nag,tend);
global abz;                   abz                      = zeros(nag,ncoh);
global abv0;                  abv0                     = zeros(nag,1);
global dis_totv;              dis_totv                 = zeros(nag,tend);
global dis_totz;              dis_totz                 = zeros(nag,ncoh);
global dis_totv0;             dis_totv0                = zeros(nag,1);
global cGv;                   cGv                      = zeros(nag,tend);
global cGz;                   cGz                      = zeros(nag,ncoh);
global cGv0;                  cGv0                     = zeros(nag,1);
global ellv;                  ellv                     = zeros(nag,tend);
global ellz;                  ellz                     = zeros(nag,ncoh);
global ellv0;                 ellv0                    = zeros(nag,1);
global gamv;                  gamv                     = zeros(nag,tend);
global gamz;                  gamz                     = zeros(nag,ncoh);
global gamv0;                 gamv0                    = zeros(nag,1);
global ivv;                   ivv                      = zeros(nag,tend);
global ivz;                   ivz                      = zeros(nag,ncoh);
global ivv0;                  ivv0                     = zeros(nag,1);
global lambdav;               lambdav                  = zeros(nag,tend);
global lambdaz;               lambdaz                  = zeros(nag,ncoh);
global lambdav0;              lambdav0                 = zeros(nag,1);
global notretv;               notretv                  = zeros(nag,tend);
global notretz;               notretz                  = zeros(nag,ncoh);
global notretv0;              notretv0                 = zeros(nag,1);
global pv;                    pv                       = zeros(nag,tend);
global pz;                    pz                       = zeros(nag,ncoh);
global pv0;                   pv0                      = zeros(nag,1);
global pcv;                   pcv                      = zeros(nag,tend);
global pcz;                   pcz                      = zeros(nag,ncoh);
global pcv0;                  pcv0                     = zeros(nag,1);
global rv;                    rv                       = zeros(nag,tend);
global rz;                    rz                       = zeros(nag,ncoh);
global rv0;                   rv0                      = zeros(nag,1);
global tauCv;                 tauCv                    = zeros(nag,tend);
global tauCz;                 tauCz                    = zeros(nag,ncoh);
global tauCv0;                tauCv0                   = zeros(nag,1);
global tauWv;                 tauWv                    = zeros(nag,tend);
global tauWz;                 tauWz                    = zeros(nag,ncoh);
global tauWv0;                tauWv0                   = zeros(nag,1);
global taulv;                 taulv                    = zeros(nag,tend);
global taulz;                 taulz                    = zeros(nag,ncoh);
global taulv0;                taulv0                   = zeros(nag,1);
global thetav;                thetav                   = zeros(nag,tend);
global thetaz;                thetaz                   = zeros(nag,ncoh);
global thetav0;               thetav0                  = zeros(nag,1);
global wv;                    wv                       = zeros(nag,tend);
global wz;                    wz                       = zeros(nag,ncoh);
global wv0;                   wv0                      = zeros(nag,1);
global yv;                    yv                       = zeros(nag,tend);
global yz;                    yz                       = zeros(nag,ncoh);
global yv0;                   yv0                      = zeros(nag,1);