classdef firm
    
    methods(Static)

        % production function and marginal products
        function out    = fY(K,L,TFP)
            global alpha;
            out = TFP.*K.^alpha.*L.^(1-alpha);
        end
        
        function out    = MPK(K,L,TFP)
            global alpha;
            out = TFP.*alpha.*(K./L).^(alpha-1);
        end
        
        function out    = MPKinv(MPKin,L,TFP)
            global alpha;
            out = L.*(MPKin./(TFP.*alpha)).^(1/(alpha-1));
        end

        function out    = MPL(K,L,TFP)
            global alpha;
            out = TFP.*(1-alpha).*(K./L).^alpha;
        end

        % invert the firm problem: finds r for given V (and LS, TFP and taxes)
        function r2 = rdemand(assetsupply, maxiter, tol, verbose)

            global K uck r qTob LD TFP tauprof delta tend;

            K2      = K;
            uck2    = uck;
            r2      = r;
            qTob2   = qTob;
            
            error   = Inf;
            iter    = 0;
            
            while true
            
                iter           = iter + 1;
                error_old      = error;
                qTob_old       = qTob2;
                
                K2             = assetsupply./qTob2;
                %K2(1)          = K(1); %predetermined
                uck2           = firm.MPK(K2,LD,TFP);
                qTob2          = (1-tauprof).*uck2 + tauprof*delta + (1-delta);
                
                error = sum(abs(qTob2-qTob_old));
                
                if (verbose)
                    fprintf("Iteration:\t%3u\t\tError:\t%8.5f\n",iter,error);
                end
                
                if (iter > maxiter)    
                    if (verbose) 
                        fprintf("No convergence!!\n"); 
                    end
                    break;
                end
                if (error < tol)       
                    if (verbose) 
                        fprintf("Convergence!!\n");
                    end
                    break; 
                end
                if (error > error_old) 
                    if (verbose) 
                        fprintf("Increasing error: stop at previous step!\n"); 
                    end
                    qTob2 = qTob_old;
                    break;
                end
            
            end
            
            r2(1:(tend-1)) = qTob2(2:tend)-1;
            r2(tend) = r2(tend-1);

        end

    end

end
