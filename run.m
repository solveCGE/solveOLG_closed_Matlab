%-------------------------------------------------------------------------%
% "solveOLG_closed_Matlab"
% 
% Solves an AK-OLG-Model, closed economy, with income effects in Matlab
% Philip Schuster, July, 2023
%
% run.m: main script, shock definition, then computes transition path
%-------------------------------------------------------------------------%

clearvars; clc; close all;

fprintf("\nSIMPLE AUERBACH-KOTLIKOFF CLOSED ECONOMY MODEL IN MATLAB\n\n");

global tend nag budget_bal;

% control center
tend            = 300;   % number of periods
nag             = 100;   % number of age groups (nag = x => max age = x-1)
budget_bal      = 3;     % budget closing instrument (1.. tauW, 2.. tauF, 3.. tauC, 4.. taul, 5.. tauprof, 6.. cG) to fulfill given debt path

% initializes all variables globally
initdata;

% run calibration routine
calib;    % <- change parameters in calib.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  POLICY SHOCK SECTION %%
%==========================
% Note: it is typically easiest to introduce shocks to period-view variables 
%       and then convert them to cohort-view variables using fun.per2coh()

%% tax shocks
%tauprof = tauprof*0+0.15;                      % profit tax rate is set to 15%
%tauprof(10:tend) = tauprof(10:tend)*0+0.15;    % profit tax rate is set to 15%, starting in 10 years (announced today)
%tauprof(1:10) = tauprof(1:10)*0+0.15;          % profit tax rate is set to 15% temporarily for next 10 years
%tauWv = tauWv*1.02; tauWz = fun.per2coh(tauWv); % wage tax is increased by 2%

%% delayed retirement (uncomment whole block)
% rag(1:10) = rag0:2/9:rag0+2; rag(11:tend)=rag0+2; 
%   notretv(:,:) = 0;
%   for tt = 1:tend
%       notretv(1:floor(rag(tt)),tt) = 1;
%       notretv(floor(rag(tt))+1,tt) = rag(tt)-floor(rag(tt));
%   end
%   notretz = fun.per2coh(notretv);                    % effective retirement age is increased linearly by 2 years over the next 10 years

%% pension cut
%pv = pv*0.95; pz = fun.per2coh(pv); % pensions are cut by 5%

%% productivity shocks
%thetav = thetav*1.02; thetaz = fun.per2coh(thetav);                       % individual productivity increases permanently by 2%
%thetav(:,1:30) = thetav(:,1:30)*1.02; thetaz = fun.per2coh(thetav);       % individual productivity increases by 2% in the next 30 years
%TFP = TFP*1.02;                                                           % total factor productivity increases permanently by 2%

%% fertility shocks
%NB = NB*1.02;            % 2% more newborns every year
NB(1:30) = NB(1:30)*1.02; % 2% more newborns every year over next 50 years

%% mortality shocks
gamv(60:nag,:) = 1-(1-gamv(60:nag,:))*0.9; gamv(nag,:)=0; gamz = fun.per2coh(gamv);  % reduction of old-age mortality by 10%

%% shock to the initial capital stock
%K(1) = K0*0.99;                                                       % 1% of capital stock is lost in first period

%% change in targeted debt path
%DG(5:20) = DG0:DG0*(0.9-1)/(20-5):DG0*0.9; DG(21:tend) = DG0*0.9; % reduce public debt by 10% over 15 years starting in 5 years

%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve transition path to new steady state
algo.solveOLG(1, 200, 1e-4, [1.0,1.0,1.0,0.5,0.7]); % solveOLG(starttime, maxiter, tol, damping_factors)

% some transition plots
figure();
plot(0:tend,[r0,r]); xlabel("time"); ylabel("real interest rate");

figure();
plot(0:tend,[w0,w]); xlabel("time"); ylabel("wage rate");

y_min = min([N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0]);
y_max = max([N/N0,Y/Y0,Inv/Inv0,Cons/Cons0,CG/CG0,A/A0,P/P0]);

figure();
plot(0:tend,[1,N/N0]*100,"DisplayName","population"); xlabel("time"); ylabel("period 0 = 100"); ylim([y_min,y_max]*100); hold on;
plot(0:tend,[1,Y/Y0]*100,"DisplayName","GDP"); hold on;
plot(0:tend,[1,Inv/Inv0]*100,"DisplayName","investment"); hold on;
plot(0:tend,[1,Cons/Cons0]*100,"DisplayName","consumption"); hold on;
plot(0:tend,[1,CG/CG0]*100,"DisplayName","public consumption"); hold on;
plot(0:tend,[1,A/A0]*100,"DisplayName","aggregate assets"); hold on;
plot(0:tend,[1,P/P0]*100,"DisplayName","pension expend.");
legend show


