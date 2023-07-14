classdef algo
    
    methods(Static)

        % computes demographic transition (and updates intervivo transfers accordingly)
        function [] = compdemo()
            
            global ivz Nv NB nag gamv Nz N Nc ivv tend fag

            % compute demography transition
            for tt = 2:tend
                Nv(1,tt)  = NB(tt);
                for i = 2:nag
                  Nv(i,tt) = Nv(i-1,tt-1)*gamv(i-1,tt-1);
                end
            end
            
            Nz = fun.per2coh(Nv);
            N  = fun.aggcoh2per(Nz);
            Nc = sum(Nv(1:(fag-1),:),1);
            
            % Compute neutral intervivo-transfers by rescaling received transfers
            for tt = 1:tend
                ivgiven          = -sum(Nv(:,tt).*ivv(:,tt).*(ivv(:,tt)<0));
                ivreceived       = sum(Nv(:,tt).*ivv(:,tt).*(ivv(:,tt)>0));
                
                ivv(ivv(:,tt)>0,tt) = ivv(ivv(:,tt)>0,tt)*(ivgiven/ivreceived);
                
                if (abs(sum(ivv(:,tt).*Nv(:,tt)))>1e-10)
                    error("ERROR IN RECOMPDEMO: Unbalanced intervivo transfers!");
                end
            
            end
            
            ivz = fun.per2coh(ivv);

        end

        % computes the further life expectancy
        function lifeexpectageN = lifeexpect(gamv)
          
          s   = size(gamv);
          nag = s(1);
          
          lifeexpectageN = zeros(nag,1);
          for a = (nag-1):-1:1
            lifeexpectageN(a) = (lifeexpectageN(a+1)+1)*gamv(a)+gamv(a+1)*(1-gamv(a));
          end
        
        end

        % main routine that solves the transition path of the full model
        function [] = solveOLG(starttime, maxiter, tol, damping_factors)
          
            global tend uck r delta tauprof K LD Inv TFP Y qTob nag w wv wz V Div TaxF tauF

            global HH_nonconvz;

            % recommended: [1.0,1.0,1.0,0.5,0.7]
            damping_budget     = damping_factors(1);
            damping_assets     = damping_factors(2);
            damping_ab         = damping_factors(3);
            damping_r          = damping_factors(4);
            damping_new_assets = damping_factors(5);
            
            fprintf("\nRunning Tatonnement Algorithm for Transition:\n\n");
            
            tic_total = tic;
            
            scaleA            = 1.0; % initialize
            scaleab           = 1.0; % initialize
            
            %===== demography ======%
            algo.compdemo(); % recomputes demographic transition
            
            for iter = 1:maxiter
                tic_iter = tic;
                
                %===== solve the firm problem for given labor demand ======%
                uck((starttime+1):tend)    = (r(starttime:(tend-1))+delta*(1-tauprof((starttime+1):tend)))./(1-tauprof((starttime+1):tend));
                K((starttime+1):tend)      = firm.MPKinv(uck((starttime+1):tend),LD((starttime+1):tend),TFP((starttime+1):tend));
                Inv(starttime:(tend-1))    = K((starttime+1):tend) - (1-delta)*K(starttime:(tend-1));
                
                Inv(tend)      = delta*K(tend);
                qTob           = (1-tauprof).*firm.MPK(K,LD,TFP) + tauprof*delta + (1-delta);
                
                Y              = firm.fY(K,LD,TFP);
                w              = firm.MPL(K,LD,TFP)./(1+tauF);
                wv             = kron(ones(nag,1),w);
                wz             = fun.per2coh(wv);
                V              = qTob.*K;
                TaxF           = tauprof.*(Y-(1+tauF).*w.*LD-delta*K);
                Div            = Y-(1+tauF).*w.*LD-Inv-TaxF;
                
                %===== solve the households' problem for given prices and tax rates ======%
                hh.HHall(starttime, (iter == 1), scaleA);
                
                global Cons Consz Nz LS notretz ellz thetaz A Az ab abz iv ivz Nw Nr;
                global P tauW TaxP Taxl Rev CG Exp tauWz pz taulz tauC cGv Nv PB DG;
                global edy edg edl eda ediv edab edw gamz Savz;

                %===== aggregation ======%
                Cons      = fun.aggcoh2per(Consz.*Nz);
                LS        = fun.aggcoh2per(notretz.*ellz.*thetaz.*Nz);
                A         = fun.aggcoh2per(Az.*Nz);
                ab        = fun.aggcoh2per(abz.*Nz);
                iv        = fun.aggcoh2per(ivz.*Nz); % should be 0 by construction
                Nw        = fun.aggcoh2per(notretz.*Nz);
                Nr        = fun.aggcoh2per((1-notretz).*Nz);
                
                % government budget
                P         = fun.aggcoh2per((1-notretz).*pz.*Nz);
                tauW      = fun.aggcoh2per(tauWz.*notretz.*ellz.*thetaz.*Nz)./LS;
                TaxP      = fun.aggcoh2per((1-notretz).*tauWz.*pz.*Nz);
                Taxl      = fun.aggcoh2per(taulz.*Nz);
                Rev       = TaxF+(tauF.*LD+tauW.*LS).*w+Taxl+tauC.*Cons+TaxP;
                CG        = sum(cGv.*Nv,1);
                Exp       = CG+P;
                
                % follow given debt-path
                PB(starttime:(tend-1))  = DG(starttime:(tend-1))-DG((starttime+1):tend)./(1+r(starttime:(tend-1)));
                PB(tend)                = r(tend)*DG(tend)/(1+r(tend));
                
                %===== excess demands ======# 
                edy       = Inv+Cons+CG-Y;
                edg       = Rev-Exp-PB;
                edl       = LD-LS;
                eda       = DG+V-A;
                ediv      = -iv;
                edab      = fun.aggcoh2per((1-gamz).*Savz.*Nz)-ab;
                edw       = 1*edy + w.*edl + ediv + edab + edg + eda - [eda(2:tend),eda(tend)]./(1+r); % Walras' Law

                % check Walras' Law: this always has to hold (even out of equilibrium)! If not there is something wrong with accounting in the model
                if (max(abs(edw(starttime:(tend-1))))> 1e-10)
                    error("Error: Walras Law does not hold!\n");
                end

                toc_iter = toc(tic_iter);
                %===== checking error and breaking loop ======% 	
                err             = sum(abs(edy(starttime:tend)))+sum(abs(edg(starttime:tend)))+sum(abs(edl(starttime:tend)))+sum(abs(eda(starttime:tend)))+sum(abs(ediv(starttime:tend)))+sum(abs(edab(starttime:tend)));
                err2            = log(err/tol);

                fprintf("Iter: %3u  scaleA: %.4f  scaleab: %.4f  non-conv.HH: %2u  Time: %5.2f sec  log of err/tol: %8.5f\n",iter, scaleA, mean(scaleab), sum(HH_nonconvz,"all"), toc_iter, err2);

                if err2 < 0.0
                  fprintf(join(repmat(' ',1,95)) + "Convergence!\n\n");
                  break;
                end
                if iter == maxiter
                  fprintf(join(repmat(' ',1,95)) + "No Convergence!\n\n");
                  break;
                end
                
                HH_nonconvz(:,:) = 0; % reset convergence counter

                %======= updating for next iteration =======%
                
                global budget_bal tauWv tauCv tauCz taul taulv fag N Nc abv rz;
                
                % budget rules
                budget_surplus  = edg*damping_budget;
                
                if budget_bal == 1
                  tauWv     = tauWv - kron(budget_surplus./(w.*LS),ones(nag,1));
                  tauWz     = fun.per2coh(tauWv);
                end
                if (budget_bal == 2)
                  tauF       = tauF - budget_surplus./(w.*LD); 
                end
                if (budget_bal == 3)
                  tauC       = tauC - budget_surplus./Cons; 
                  tauCv      = kron(tauC,ones(nag,1));
                  tauCz      = fun.per2coh(tauCv);
                end
                if (budget_bal == 4)
                  taul             = taul - budget_surplus./(N-Nc); % no updating of taulv, is this fixed now?
                  taulv(fag:nag,:) = kron(taul,ones(nag-fag+1,1));
                  taulz            = fun.per2coh(taulv);
                end
                if (budget_bal == 5)
                  tauprof    = tauprof - budget_surplus./(Y-(1+tauF).*w.*LD-delta.*K); 
                end
                if (budget_bal == 6)
                  cGv        = cGv + kron(budget_surplus./N,ones(nag,1));
                  CG         = sum(cGv.*Nv,1);
                end

                % price updating
                newassets       = damping_new_assets*(A-DG) + (1-damping_new_assets)*V;
                r_new           = firm.rdemand(newassets, 20, 1e-6, false);    % rdemand(assetsupply, maxiter, tol, verbose)
                r               = damping_r*r_new + (1-damping_r)*r;
                rv              = kron(ones(nag,1),r);
                rz              = fun.per2coh(rv);
            
                scaleab         = 1+(fun.aggcoh2per((1-gamz).*Savz.*Nz)./ab-1)*damping_ab;
                abv             = abv.*kron(scaleab,ones(nag,1)); % matrix??? dim of scaleab??
                abz             = fun.per2coh(abv);
                LD              = LS;
                scaleA          = 1+((DG(starttime)+V(starttime))/A(starttime)-1)*damping_assets;
                
            end

            global Av Consv lambdav lambdaz Savv dis_totv dis_totz ellv pcv pcz yv yz;

            % convert cohort-view variables back to period-view variables
            % (those where only cohort-view variables were altered in solveOLG)
            Av       = fun.coh2per(Az);
            Consv    = fun.coh2per(Consz);
            lambdav  = fun.coh2per(lambdaz);
            Savv     = fun.coh2per(Savz);
            dis_totv = fun.coh2per(dis_totz);
            ellv     = fun.coh2per(ellz);
            pcv      = fun.coh2per(pcz);
            yv       = fun.coh2per(yz);
            
            toc_total = toc(tic_total);
            fprintf("Computation time:\t%.3f sec\n", toc_total);
            fprintf("CHECK SOLUTION:\t\t%.5f\n",sum(abs(edy)+abs(edl)+abs(edg)+abs(eda)+abs(ediv)+abs(edab)));

        end



    end
end
