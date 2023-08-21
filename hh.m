classdef hh
    methods(Static)
        
        % solves the household problem for cohort z
        function out = HH_root(lambdain, sage, z) 

            global lambdaz rho gamz rz pcz tauCz Consz fag nag sigma ellz;
            global wz tauWz thetaz parlv0 sigL dis_totz parlv1 yz notretz pz taulz;
            global Az ivz abz Savz;
            
            % EULER EQUATION: solve forward in age
            lambdaz(sage,z)    = lambdain;
            if sage < nag
                for a = sage:(nag-1)
                    lambdaz(a+1,z) = lambdaz(a,z)/((1/(1+rho))*gamz(a,z)*(1+rz(a,z)));
                end
            end
            
            % CONSUMPTION
            pcz(sage:nag,z)      = 1+tauCz(sage:nag,z);
            Consz(sage:nag,z)     = (pcz(sage:nag,z).*lambdaz(sage:nag,z)).^(-sigma);
            
            % HOURS SUPPLY
            ellz(sage:nag,z)     = ((wz(sage:nag,z).*(1-tauWz(sage:nag,z)).*thetaz(sage:nag,z)./pcz(sage:nag,z).*(Consz(sage:nag,z).^(-1/sigma)))./parlv0(sage:nag)).^sigL;
            dis_totz(sage:nag,z) = (sigL/(1+sigL))*parlv0(sage:nag).*ellz(sage:nag,z).^((1+sigL)/sigL)-parlv1(sage:nag);
            
            % CONSUMPTION AND SAVINGS
            yz(sage:nag,z)       = notretz(sage:nag,z).*(wz(sage:nag,z).*(1-tauWz(sage:nag,z)).*ellz(sage:nag,z).*thetaz(sage:nag,z))+(1-notretz(sage:nag,z)).*(1-tauWz(sage:nag,z)).*pz(sage:nag,z)-taulz(sage:nag,z);
            
            % ASSETS: solve forward in age
            Az(1,z)         = 0;
            
            if sage < nag
                for a = sage:(nag-1)
                    Az(a+1,z)   = (1+rz(a,z)).*(Az(a,z)+yz(a,z)+ivz(a,z)+abz(a,z)-pcz(a,z).*Consz(a,z)); % if sage > 1 take previous age entry in Az as starting value! (i.e. has to be given globally not passed in function)
                end
            end
            Savz(sage:nag,z)  = Az(sage:nag,z)+yz(sage:nag,z)+ivz(sage:nag,z)+abz(sage:nag,z)-pcz(sage:nag,z).*Consz(sage:nag,z);
            
            out = Savz(nag,z);
        end

        % solves the household problem for cohort z given an initial guess for lamdba
        function [] = HH(sage, z, maxiter, stol, atol)
            
            global lambdaz HH_nonconvz nag Savz;

            err            = Inf;
            iter           = 0;
            trys           = 0;
            stepsize       = 1e-6; % for numerical gradient
            
            lambdatrys     = [1.0,0.5,1.5,0.25,1.25,0.1,1.0];
            s              = size(lambdatrys);
            maxtrys        = s(2);
            while_continue = true;
            
            while (while_continue)
            
                while_continue = false;
                lambdazsave    = lambdaz(sage,z);
                
                while (((err > stol)||(abs(Savz(nag,z)) > atol)) && (trys < maxtrys))
                  
                    trys = trys + 1;
                    iterpertry = 0;
                    lambdaz1 = lambdazsave*lambdatrys(trys);
                    
                    breakwhile = false;
                    while ((err > stol) && (iterpertry < maxiter) && (~breakwhile))
                        if (iterpertry == 0) % Newton step for first iteration
                            f2 = hh.HH_root(lambdaz1+stepsize,sage,z); iter = iter + 1;
                            
                            if (~isfinite(f2))
                                breakwhile = true; break;
                            end
                            f1 = hh.HH_root(lambdaz1,sage,z); iter = iter + 1;
                            
                            if (~isfinite(f1))
                                breakwhile = true; break;
                            end
                            lambdaz2 = lambdaz1 - f1*stepsize/(f2-f1);
                            if (~isfinite(lambdaz2)||(lambdaz2<0))
                                breakwhile = true; break;
                            end
                        else % Secant method
                            f1 = hh.HH_root(lambdaz1,sage,z); iter = iter + 1;
                        
                            if (~isfinite(f1))
                                breakwhile = true; break;
                            end
                            lambdaz2 = lambdaz1 - f1*(lambdaz1-lambdaz0)/(f1-f0);
                            if (~isfinite(lambdaz2)||(lambdaz2<0))
                                breakwhile = true; break;
                            end
                        end
                        err = abs(lambdaz2-lambdaz1);
                        lambdaz0 = lambdaz1;
                        lambdaz1 = lambdaz2;
                        f0       = f1;
                        iterpertry = iterpertry + 1;
                    end
    
                end
            end
            
            if (abs(Savz(nag,z)) > atol)
                HH_nonconvz(nag,z) = 1; % counter
            end

        end

        % solves the household problem for all cohorts
        function [] = HHall(starttime, calibinit, scaleA)

            global ncoh nag fag Az Av0;

            for z = starttime:ncoh
                if z <= (nag-fag+starttime-1)
                    if calibinit
                        Az(:,z) = Av0;
                    end
                    Az(nag-(z-starttime),z) = Az(nag-(z-starttime),z)*scaleA;
                    hh.HH(nag-(z-starttime), z, 30, 1e-10, 0.1); % HH(sage, z, maxiter, stol, atol)
                else
                    hh.HH(fag, z, 30, 1e-10, 0.1);
                end
            end
        end
    
    end
end

