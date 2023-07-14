
classdef fun
    
    methods(Static)

        % formatted reporting
        function [] = report(reporttext,reportcalc)
          cursorstart     = 45;
          %cursorstart     = max(strlength(reporttext),cursorstart)+3*(reportcalc>=0)+2*(reportcalc<0);
          cursorstart     = max(strlength(reporttext),cursorstart)+3;
          charfill        = join(repmat(' ',1,cursorstart-strlength(reporttext)));
          
          fprintf(reporttext + charfill + '%12.8f\n',reportcalc);
          
        end

        % convert variable from cohort-view to period-view
        function outmat = coh2per(inmat)
          
            s = size(inmat);
            maxage = s(1);
            numcoh = s(2);
            numper = numcoh-(maxage-1);
            
            if (numper<=0), error("coh2per: insufficient number of columns in input matrix"); end
            
            outmat = zeros(maxage, numper);
        
            for a = 1:maxage
                outmat(a,:) = inmat(a,(maxage:numcoh)-(a-1));
            end
        end
        
        % convert variable from period-view to cohort-view
        function outmat = per2coh(inmat)
            
            s = size(inmat);
            maxage = s(1);
            numper = s(2);
            numcoh = numper+(maxage-1);
        
            calibvec = inmat(:,1);
            
            outmat = zeros(maxage,numcoh);
            
            for a = 1:maxage
                if (a < maxage), outmat(a,1:(maxage-a)) = calibvec(a)*ones(maxage-a,1); end
                outmat(a,(maxage:numcoh)-(a-1)) = inmat(a,:);
                if (a > 1), outmat(a,(numcoh-(a-2)):numcoh) = inmat(a,numper)*ones(a-1,1); end
            end
               
        end
        
        % convert variable from cohort-view to period-view and aggregate over age
        function outmat = aggcoh2per(inmat) 
            outmat = sum(fun.coh2per(inmat),1);
        end

        % simple Newton method
        function [xout, xeval, flag] = fNewton(func, x0, tol)
            
            maxiter = 30;
            flag    = 0;

            for i = 1:maxiter
                ff = func(x0);
                if (max(abs(ff))<tol)
                    flag = 1;
                    break;
                end
                x0 = x0 - (fun.fJac(func,x0)\ff)';
            end
            
            xout  = x0;
            xeval = func(x0);
        end

        % Jacobian matrix
        function Jac = fJac(func, x0)
        
            stepsize = 1e-6;
        
            s = size(x0);
            Jac = zeros(s(2),s(2));
        
            for i = 1:s(2)
                eps = zeros(1,s(2));
                eps(i) = stepsize;
                f2 = func(x0+eps);
                f1 = func(x0);
                Jac(:,i) = (f2-f1)./stepsize;
            end

        end

    end

end