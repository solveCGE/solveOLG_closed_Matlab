function retvar = calibfind(xcalib0)
  
    global nag;
    global fag;                   
    global rho;                   
    global sigma;

    global A0;                    
    global Cons0;                 
    global CG0;                   
    global DG0;                   
    global Exp0;                  
    global Inv0;                  
    global LD0;                   
    global LS0;                   
    global N0;                    
    global Nc0;                   
    global Nr0;                   
    global Nw0;                   
    global P0;                    
    global PB0;                   
    global Rev0;                  
    global TaxF0;                 
    global V0;                    
    global Y0;                    
    global ab0;                   
    global eda0;                  
    global edab0;                 
    global edg0;                  
    global ediv0;                 
    global edl0;                  
    global edy0;                  
    global iv0;                   
    global r0;
    global tauC0;                 
    global tauF0;                 
    global tauW0;                 
    global taul0;                 
    global w0;

    global Av0;                   
    global Consv0;                
    global Nv0;                   
    global Savv0;                 
    global abv0;                  
    global cGv0;                  
    global ellv0;                 
    global gamv0;                 
    global ivv0;                  
    global lambdav0;              
    global notretv0;              
    global pv0;                   
    global pcv0;                  
    global rv0;                   
    global tauWv0;                
    global taulv0;                
    global thetav0;               
    global wv0;                   
    global yv0;

    global cGv0_profile;
    
    retvar = zeros(5,1);
    
    rho       = xcalib0(1);
    cGscale   = xcalib0(2);
    taul0     = xcalib0(3);
    ab0       = xcalib0(4);
    lambdain  = xcalib0(5);
    
    abv0(fag:nag)    = ab0/(N0-Nc0)*ones(nag-fag+1,1); % children do not receive accidental bequest (workers start out with 0 assets)
    taulv0(fag:nag)  = taul0;
    cGv0             = cGv0_profile+cGscale;
    
    % INCOME
    yv0     = notretv0.*(wv0.*(1-tauWv0).*ellv0.*thetav0)+(1-notretv0).*(1-tauWv0).*pv0-taulv0;
    
    % CONSUMPTION FOR ALL AGE GROUPS
    
    % Euler equation
    lambdav0(fag) = lambdain;
    for a = fag:(nag-1)
        lambdav0(a+1) = lambdav0(a)/((1/(1+rho))*gamv0(a)*(1+rv0(a)));
    end
    Consv0(fag:nag) = (pcv0(fag:nag).*lambdav0(fag:nag)).^(-sigma);
    
    % assets
    Av0(fag) = 0;
    for a = (fag+1):nag
     Av0(a)     = (1+rv0(a-1))*(Av0(a-1)+yv0(a-1)+ivv0(a-1)+abv0(a-1)-pcv0(a-1)*Consv0(a-1));
    end
    Savv0   = Av0+yv0+ivv0+abv0-pcv0.*Consv0;
    
    % AGGREGATION
    A0      = sum(Av0.*Nv0);                                % total assets
    P0      = sum((1-notretv0).*pv0.*Nv0);                  % expend pensions
    CG0     = sum(cGv0.*Nv0);                               % government consumption
    Exp0    = CG0+P0;                                       % total primary expenditures
    tauW0   = sum(tauWv0.*notretv0.*ellv0.*thetav0.*Nv0)/LS0;  % average wage tax rate
    Rev0    = TaxF0+(tauF0*LD0+tauW0*LS0)*w0+taul0*(Nw0+Nr0)+tauC0*Cons0+sum((1-notretv0).*tauWv0.*pv0.*Nv0); % total revenues
    PB0     = DG0*r0/(1+r0);                                % primary balance
    
    % EXCESS DEMANDS
    edy0    = Cons0+CG0+Inv0-Y0;      % goods market
    edl0    = LD0-LS0;                % labor market
    eda0    = DG0+V0-A0;              % assets market
    edg0    = Rev0-Exp0-PB0;          % government budget
    ediv0   = -iv0;                   % intervivo transfers resource constraint
    edab0   = sum((1-gamv0).*Savv0.*Nv0)-ab0; % accidental bequest resource constraint
    
    retvar(1)   = edy0;
    retvar(2)   = edg0;
    retvar(3)   = sum(Consv0.*Nv0)-Cons0;
    retvar(4)   = edab0;
    retvar(5)   = Savv0(nag);
    
end
