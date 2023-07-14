fprintf("Start calibration: \n");

% parameters
delta                    = 0.05;                 % depreciation rate
r0                       = 0.04;                 % real interest rate
sigma                    = 0.9;                  % elasticity of inter-temporal substitution
sigL                     = 0.3;                  % labor supply elasticity

% note: ages are off-set by 1 year, e.g. age group 1 contains 0-year olds
fag                      = 14;                   % first economically active age-group (age 15)
rag0                     = 61.3;                 % retirement age group (retirement age 62.3), non-whole numbers allowed
iag0                     = 51;                   % first age group giving inter-vivo transfers

ivpc                     = 0.2;                  % intervivo transfer received per capita

% some normalizations
N0                       = 100.0;                % population
Y0                       = 100.0;                % GDP
L0                       = 30.0;                 % total labor supply in efficiency units
w0                       = 2.0;                  % wage rate
Cons0                    = 0.55*Y0;              % consumption share (calibrated using taul0)

%%% DEMOGRAPHY %%%
for i = 1:nag
  gamv0(i) = 1-0.89^(nag-i+1);                   % some simple profile
end
  
% survival of last age group is 0
gamv0(nag) = 0;
  
% compute demography
Nv0(1)  = 1;
for i = 2:nag
  Nv0(i) = Nv0(i-1)*gamv0(i-1);
end
  
% rescale population
NB0             = 1/sum(Nv0)*N0;
Nv0             = Nv0/sum(Nv0)*N0;

avage0 = sum(Nv0.*(0:(nag-1))')/N0;
fun.report("REPORT: Average age:",avage0);
lifeexpect0 = algo.lifeexpect(gamv0);
fun.report("REPORT: Life-expectancy at age 0:", lifeexpect0(1));
fun.report("REPORT: Life-expectancy at age 65:", lifeexpect0(66));

%%% AGE PROFILES %%%

% indicator for not-retired
notretv0(1:floor(rag0))   = 1;                       % not retired
notretv0(floor(rag0)+1)   = rag0-floor(rag0);        % partly retired

% intervivo-transfers
%ivv0(iag0:nag)            = -seq(from=ivpc,to=ivpc*2,length.out=nag-iag0+1); % some increasing profile (from ivpc to 2*ivpc)
ivv0(iag0:nag)            = -(ivpc:(ivpc*2-ivpc)/(nag-iag0):ivpc*2); % some increasing profile (from ivpc to 2*ivpc)
ivv0(fag:(iag0-1))        = -sum(ivv0(iag0:nag).*Nv0(iag0:nag))/sum(Nv0(fag:(iag0-1)))*ones(iag0-fag,1);

iv0                       = sum(ivv0.*Nv0);
if (abs(iv0)>1e-10), error("ERROR: UNBALANCED INTERVIVO TRANSFERS!"); end

thetav0                     = zeros(nag,1);                                % labor productivity parameters
theta_peak                  = floor(rag0)-10;                              % assumption: productivity peaks 10 years before retirement
thetav0(fag:theta_peak)     = 0.7:(1-0.7)/(theta_peak-fag):1;
thetav0((theta_peak+1):nag) = 1:(0.1-1)/(nag-theta_peak-1):0.1;

ellv0                       = L0/sum(Nv0.*thetav0.*notretv0)*ones(nag,1);   % labor supply

% partition of population
Nc0   = sum(Nv0(1:(fag-1)));        % number of children
Nw0 	= sum(notretv0.*Nv0)-Nc0; 		  % number of workers
Nr0 	= sum((1-notretv0).*Nv0); 	    % number of retirees
fun.report("REPORT: Old-age dependency ratio:",sum(Nv0(66:nag))/sum(Nv0(16:65)));
fun.report("REPORT: Economic dependency ratio:",(Nc0+Nr0)/Nw0);
fun.report("CHECK: Newborns - deaths:", sum((1-gamv0).*Nv0)-NB0);
fun.report("CHECK: Children + workers + retriees - pop.:", Nc0+Nw0+Nr0-N0);

%%% POLICY PARAMETERS %%%
tauWv0                   = 0.15*ones(nag,1);             % wage tax rate worker & retiree
tauF0                    = 0.2;                          % payroll tax rate
tauC0                    = 0.2;                          % consumption tax rate
tauprof0                 = 0.1;                          % profit tax rate
pv0                      = 0.65*sum(w0.*ellv0.*thetav0.*Nv0)/N0*ones(nag,1);  % old-age pension (65% of average wage earnings)
DG0                      = 0.6*Y0;                       % government debt level (60% of GDP)

% cGv0 is used to balance budget in calibration
global cGv0_profile;
cGv0_profile             = 0.2*ones(nag,1);
cGv0_profile(1:25)       = 0.4:(0.2-0.4)/24:0.2;
cGv0_profile(55:nag)     = 0.2:(1.0-0.2)/(nag-55):1.0; % some U-shaped profile

% price of consumption and age specific prices and tax rates (but the same for all age groups)
pc0     = 1+tauC0;
tauCv0  = tauC0*ones(nag,1);
pcv0    = pc0*ones(nag,1);
wv0     = w0*ones(nag,1);
rv0     = r0*ones(nag,1);

LS0     = sum(notretv0.*ellv0.*thetav0.*Nv0); % aggregate labor supply
LD0     = LS0;
uck0    = (r0+delta*(1-tauprof0))/(1-tauprof0); % user-cost of capital
K0      = (Y0-(1+tauF0)*w0*LD0)/uck0;
Inv0    = delta*K0;
alpha   = K0*uck0/(K0*uck0+LS0*((1+tauF0)*w0));
qTob0   = (1-tauprof0)*alpha*Y0/K0 + tauprof0*delta + (1-delta); % = 1+r0
TFP0    = Y0/((K0^alpha)*(LS0^(1-alpha)));
%LD0     = ((1-alpha)*TFP0/((1+tauF0)*w0))^(1/alpha)*K0; % also true
TaxF0   = tauprof0*(Y0-(1+tauF0)*w0*LD0-(delta*K0));
Div0    = Y0-(1+tauF0)*w0*LD0-Inv0-TaxF0;
V0      = (1+r0)*Div0/r0;

% MATCH CALIBRATION TARGETS;
xcalib0 = [0.01, 0.3719, 0.40, 13, 1];

[xout, xeval, exitflag] = fun.fNewton(@calibfind, xcalib0, 1e-10);
if exitflag~= 1, error('Calibration could not be found!\n'); end

%%% CALIBRATION OF LABOR SUPPLY MARGINS %%%

% set parl0 in order to reproduce ell0, FOC ell0
parlv0  = (wv0.*(1-tauWv0).*thetav0./pcv0).*(ellv0.^(-1/sigL)).*(Consv0.^(-1/sigma)); parlv0(1:(fag-1)) = 0;
% set parl1 in order to normalize disutility of labor to 0
parlv1  = (sigL/(1+sigL))*parlv0.*(ellv0.^((1+sigL)/sigL));
dis_totv0= (sigL/(1+sigL))*parlv0.*(ellv0.^((1+sigL)/sigL))-parlv1;

fun.report("REPORT: Asset-to-output ratio:", A0/Y0);

checkA0         = sum(Av0.*Nv0)-A0;
checkAv0        = Av0(nag)+yv0(nag)+ivv0(nag)+abv0(nag)-pc0*Consv0(nag); % end of period assets of last age group are zero
checkN0         = sum(Nv0)-N0;

chkcalib = [edy0,edl0,edg0,ediv0,eda0,edab0,checkA0,checkAv0,checkN0];

fun.report("CHECK: Calibration:",sum(chkcalib));

% fill time-dependent variables with calibration values
Cons                          = Cons0*ones(1,tend);
DG                            = DG0*ones(1,tend);
Inv                           = Inv0*ones(1,tend);
LD                            = LD0*ones(1,tend);
LS                            = LS0*ones(1,tend);
K                             = K0*ones(1,tend);
N                             = N0*ones(1,tend);
NB                            = NB0*ones(1,tend);
PB                            = PB0*ones(1,tend);
TFP                           = TFP0*ones(1,tend);
ab                            = ab0*ones(1,tend);
pc                            = pc0*ones(1,tend);
r                             = r0*ones(1,tend);
rag                           = rag0*ones(1,tend);
tauC                          = tauC0*ones(1,tend);
tauF                          = tauF0*ones(1,tend);
tauW                          = tauW0*ones(1,tend);
taul                          = taul0*ones(1,tend);
tauprof                       = tauprof0*ones(1,tend);
uck                           = uck0*ones(1,tend);

% fill time-dependent and age-dependent variables with calibration values
Av                            = kron(Av0, ones(1,tend));
Az                            = kron(Av0, ones(1,ncoh));
Consv                         = kron(Consv0, ones(1,tend));
Consz                         = kron(Consv0, ones(1,ncoh));
Nv                            = kron(Nv0, ones(1,tend));
Nz                            = kron(Nv0, ones(1,ncoh));
Savv                          = kron(Savv0, ones(1,tend));
Savz                          = kron(Savv0, ones(1,ncoh));
abv                           = kron(abv0, ones(1,tend));
abz                           = kron(abv0, ones(1,ncoh));
cGv                           = kron(cGv0, ones(1,tend));
cGz                           = kron(cGv0, ones(1,ncoh));
ellv                          = kron(ellv0, ones(1,tend));
ellz                          = kron(ellv0, ones(1,ncoh));
gamv                          = kron(gamv0, ones(1,tend));
gamz                          = kron(gamv0, ones(1,ncoh));
ivv                           = kron(ivv0, ones(1,tend));
ivz                           = kron(ivv0, ones(1,ncoh));
lambdav                       = kron(lambdav0, ones(1,tend));
lambdaz                       = kron(lambdav0, ones(1,ncoh));
notretv                       = kron(notretv0, ones(1,tend));
notretz                       = kron(notretv0, ones(1,ncoh));
pv                            = kron(pv0, ones(1,tend));
pz                            = kron(pv0, ones(1,ncoh));
tauCv                         = kron(tauCv0, ones(1,tend));
tauCz                         = kron(tauCv0, ones(1,ncoh));
tauWv                         = kron(tauWv0, ones(1,tend));
tauWz                         = kron(tauWv0, ones(1,ncoh));
taulv                         = kron(taulv0, ones(1,tend));
taulz                         = kron(taulv0, ones(1,ncoh));
thetav                        = kron(thetav0, ones(1,tend));
thetaz                        = kron(thetav0, ones(1,ncoh));
rv                            = kron(rv0, ones(1,tend));
rz                            = kron(rv0, ones(1,ncoh));
wv                            = kron(wv0, ones(1,tend));
wz                            = kron(wv0, ones(1,ncoh));

% some optional plots of the calibration

plot((1:nag)-1,Av0); xlabel("age"); ylabel("assets");

figure();
plot((1:nag)-1,Consv0, "DisplayName", "consumption"); xlabel("age"); hold on;
plot((1:nag)-1,notretv0.*ellv0.*thetav0.*wv0.*(1-tauWv0), "DisplayName", "net labor income");
plot((1:nag)-1,(1-notretv0).*pv0.*(1-tauWv0), "DisplayName", "net pension income"); hold on;
plot((1:nag)-1,cGv0, "DisplayName", "public consumption");
legend show
