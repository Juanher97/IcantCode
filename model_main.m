clc
clear all
tic



%%                      
                        %%% 1. Parameters and Setup %%%

beta        = 0.87;
delta       = 0.06;
gamma       = 2;
alpha       = 0.35;
wage        = 1;
za          = 1; 
epsilon     = 3;
mu          = epsilon/(epsilon-1);
D           = 1;

%overhead labor;
phil            = 0.00;    %baseline model
%phil            = 0.135;  %robustness check with overhead labor

%permanent productivity differences;
nzp = 2;
zp_sd_target = 0.43; % target standard deviation;
piL=0.80;
zLspace=linspace(0.0001,0.9999,1000000);
difference=piL*(zLspace-1).^2 + (1-piL)*((piL-piL.*zLspace)./(1-piL)).^2-zp_sd_target^2;
[value,index]=min(abs(difference));
zp_l=zLspace(index);
zp_h=(1-piL*zp_l)/(1-piL);                   
zp_grid = [zp_l ,zp_h];
zpprob = eye(nzp);   

%for extension with exogenous labor wedges;
ntau = 1; %baseline model has ntau=1, extension with exogenous labor wedges has ntau=2;
if ntau==2
    tau_grid = linspace(-0.29,0.29,ntau);
    tauprob=[0.81, 0.19; 0.19, 0.81];
elseif ntau == 1
    tau_grid = 0.00;
    tauprob = 1;
end

%for extension with umeasured capital;
nq = 1; %baseline model has nq=1, extension with unmeasured capital has nq=2;
if nq==2
    qprob = eye(nq);    
    q_grid = linspace(-0.18,0.18,nq);
elseif nq == 1
    q_grid = 0.00;
    qprob = 1;
end

%transitory log productivity shocks
nzt         = 11;
rho_zt      = 0.59;
sigma_zt    = 0.13;
mu_zt       = -(sigma_zt^2)/(2*(1+rho_zt));
[zt_grid,ztprob] = tauchen(nzt,mu_zt,rho_zt,sigma_zt,3);

%interest rate process
runexp=1;  %=1, all changes in r are unexpected (baseline), =0, AR(1) process;
nr = 6;       
r_grid  =[0.01, 0.02, 0.03, 0.04, 0.05, 0.10];
if runexp==1
    rprob   = eye(nr);      
    elseif runexp==0                               
    %initial drop unexpected and then from AR(1) process between 1994-2011;        
    rho_r     = 0.50;
    sigma_r   = 0.0086;
    mu_r      = 0.03*rho_r;
    [r_grid,rprob] = tauchen(nr-1,mu_r,rho_r,sigma_r,2.014); 
    r_grid=[r_grid,max(r_grid)];
    rprob=[rprob,zeros(nr-1,1);
        0.00 0.00 0.00 0.00 0.00 1.00;];    
end

%financial frictions parameters in collateral contraint: k'< chi0*a' + chi1*(exp(k')-1);
%baseline HeF model:
chi0=0.98;
chi1=0.047;
%standard model with homogeneous frictions (HoF):
%chi0=1.06;
%chi1=0.00;
%NoF (no financial frictions) model:
%chi0=10^10;
%chi1=0;
%adjustment costs calibrated to match K response 99-07:
%chi0=0.98;
%chi1=0.047;
%extension of HeF model recalibrated to overhead labor:
%chi0=0.98;
%chi1=0.047;
%extension of HeF model recalibrated to exogenous labor wedge shocks:
%chi0=1.01;
%chi1=0.050;
%extension of HeF model recalibrated to unmeasured capital:
%chi0=1.01;
%chi1=0.037;
%extension of HeF model with AR(1) process for r:
%chi0=1.02;
%chi1=0.042;

%adjustment cost parameter;
psi = 3.2;    %baseline model HeF;
%psi = 3.2;    %model HoF;
%psi = 3.5;    %model NoF;
%psi = 7.6;    %adjustment costs calibrated to match K response 99-07;
%psi = 3.2;    %model recalibrated to overhead labor;
%psi = 2.3;    %model recalibrated to exogenous labor wedge shocks;
%psi = 1.6;    %model recalibrated to unmeasured capital;
%psi = 3.1;    %model with AR(1) process for r;

%grids for capital and net worth;
k_l=0.01;  
k_h=6.0;
nk = 120;  
k_grid = linspace(k_l,k_h,nk);

a_l=0.01;
a_h=3.0;
na = 120;     
a_grid = linspace(a_l,a_h,na);

n_choice = nk*na; 
n_state  = nzp*nzt*nr*ntau*nq;

% Create matrices for zp, zt, r, tau, q, a', and k';
r_grid_ind  = 1:nr;
zt_grid_ind = 1:nzt;
zp_grid_ind = 1:nzp;
tau_grid_ind = 1:ntau;
q_grid_ind = 1:nq;

Q = repmat(q_grid',nzp*nzt*nr*ntau,1);
Q_ind = repmat(q_grid_ind',nzp*nzt*nr*ntau,1);

TAU  = repmat(tau_grid',nq,1);
TAU = sort(TAU);
TAU  = repmat(TAU,nzp*nzt*nr,1);
TAU_ind  = repmat(tau_grid_ind',nq,1);
TAU_ind = sort(TAU_ind);
TAU_ind  = repmat(TAU_ind,nzp*nzt*nr,1);

R  = repmat(r_grid',ntau*nq,1);
R = sort(R);
R  = repmat(R,nzp*nzt,1);
R_ind  = repmat(r_grid_ind',ntau*nq,1);
R_ind = sort(R_ind);
R_ind  = repmat(R_ind,nzp*nzt,1);

ZT  = repmat(zt_grid',nr*ntau*nq,1);
ZT  = sort(ZT);
ZT  = repmat(ZT,nzp,1);
ZT_ind  = repmat(zt_grid_ind',nr*ntau*nq,1);
ZT_ind  = sort(ZT_ind);
ZT_ind  = repmat(ZT_ind,nzp,1);


ZP  = repmat(zp_grid',nzt*nr*ntau*nq,1);
ZP  = sort(ZP);
ZP_ind  = repmat(zp_grid_ind',nzt*nr*ntau*nq,1);
ZP_ind  = sort(ZP_ind);

EXOG = [ZP,ZT,R,TAU,Q];
EXOG_ind = [ZP_ind,ZT_ind,R_ind,TAU_ind,Q_ind];

EXOG = EXOG';
EXOG_ind = EXOG_ind';

% Combine a_grid and k_grid for vectorization purposes;
A = repmat(a_grid',nk,1);  
A = sort(A);               
K = repmat(k_grid',na,1);
G = [A,K];                 

%probability of transitioning from some (zt,zp,tau,r) to some (zt',zp',tau',r');
prob = zeros(n_state,n_state);

for i_state =1:n_state
    
    izp  = EXOG_ind(1,i_state);
    izt  = EXOG_ind(2,i_state);
    ir   = EXOG_ind(3,i_state);
    itau = EXOG_ind(4,i_state);
    iq   = EXOG_ind(5,i_state);
    
    for i_state_next = 1:n_state

        izpnext     = EXOG_ind(1,i_state_next);
        iztnext     = EXOG_ind(2,i_state_next);
        irnext      = EXOG_ind(3,i_state_next);
        itaunext    = EXOG_ind(4,i_state_next);
        iqnext      = EXOG_ind(5,i_state_next);
        
        prob(i_state,i_state_next)=ztprob(izt,iztnext)*zpprob(izp,izpnext)*rprob(ir,irnext)*tauprob(itau,itaunext)*qprob(iq,iqnext);
    end 
    
end

%%                      
                %%% 2. Value Function Iteration %%%

iter_vfi     = 0;
iter_howard  = 0;
error_vfi    = 10^15; 
error_vfi_old= error_vfi;
VFIcontinue  = 1;
error_tolerance = 10^(-6);
if gamma<2.5 && runexp==1
penalty     = -10^(9);
else
penalty     = -10^(12);
end
Nh          = 20; % put -1 to shut-off howard improvement;

indicator    = ones(n_choice,n_state);

pol_ind_kp   = zeros(n_choice,n_state);
pol_ind_ap   = zeros(n_choice,n_state);
V_new        = zeros(n_choice,n_state);
V_old        = penalty*ones(n_choice,n_state);
Vh_old       = penalty*ones(n_choice,n_state);
Vh_new       = penalty*ones(n_choice,n_state);

% The following lines calculate auxillary variables l, y, T, and pi, together with consumption "c" and collateral value "collateral" for each state (ia,ik,izp,izt,ir,itau); 
% Later we create an "indicator" that takes the value of one if the choice does not violate either c > 0 or collateral >= 0) for a given state;

collateral = chi0 .*G(:,1) + chi1.*(exp(G(:,2))-1)  - G(:,2); 
    
 for i_state = 1:n_state
     
    izp  = EXOG_ind(1,i_state);
    izt  = EXOG_ind(2,i_state);
    ir   = EXOG_ind(3,i_state);
    itau = EXOG_ind(4,i_state);
    iq   = EXOG_ind(5,i_state);
    
    for iak_state = 1: n_choice
        
        ik = rem(iak_state-1,nk)+1; 
        ia = floor((iak_state- 1)/nk)+1;
                        
        dummy = ones(na*nk,1); 

        l  = phil + ((za*exp(zt_grid(izt))*zp_grid(izp))^((epsilon-1)/(1+alpha*(epsilon-1))))*((k_grid(ik)+q_grid(iq))^(alpha*(epsilon-1)/(1+alpha*(epsilon-1))))*(D^(epsilon/(1+alpha*(epsilon-1))))*(mu^(-epsilon/(1+alpha*(epsilon-1))))*((wage*(1+tau_grid(itau))/(1-alpha))^(-epsilon/(1+alpha*(epsilon-1))));
        y  = (za*exp(zt_grid(izt))*zp_grid(izp))*((k_grid(ik)+q_grid(iq))^(alpha))*((l-phil)^(1-alpha));
        T  = wage*tau_grid(itau)*l+wage*(1+tau_grid(itau))*phil-(r_grid(ir) + delta )*q_grid(iq);
        pi = D*(y^((epsilon-1)/epsilon)) - wage*(1+tau_grid(itau))*l;

        c = pi - (r_grid(ir) + delta )*k_grid(ik) + (1+r_grid(ir))*a_grid(ia) + T - G(:,1) - (psi*(G(:,2) - k_grid(ik)).^2)/(2*k_grid(ik));
       
        CAPeffective = k_grid(ik)+q_grid(iq);
        
        dummy(c<=0 | collateral<0 | CAPeffective<=0 ) = -1;
        if max(dummy) < 0
            indicator(iak_state,i_state) = 0;
        end
    end
 end
                
ZP_mat = repmat(EXOG(1,:),na*nk,1);
ZT_mat = repmat(EXOG(2,:),na*nk,1);
R_mat = repmat(EXOG(3,:),na*nk,1);
Q_mat = repmat(EXOG(5,:),na*nk,1);

AP_mat = repmat(G(:,1),1,n_state);
KP_mat = repmat(G(:,2),1,n_state);

l_vfi  = zeros(n_choice,n_state);
y_vfi  = zeros(n_choice,n_state);
T_vfi  = zeros(n_choice,n_state);
pi_vfi = zeros(n_choice,n_state);

location = zeros(n_choice, n_state);

EXP      = zeros(n_choice, n_state);

for iak_state = 1:n_choice

    [ik,ia] = ind2sub([nk,na],iak_state);
    
    l_vfi(iak_state,:)  = phil + ((za*exp(EXOG(2,:)).*EXOG(1,:)).^((epsilon-1)/(1+alpha*(epsilon-1)))).*((k_grid(ik)+EXOG(5,:)).^(alpha*(epsilon-1)/(1+alpha*(epsilon-1)))).*(D^(epsilon/(1+alpha*(epsilon-1)))).*(mu^(-epsilon/(1+alpha*(epsilon-1)))).*((wage.*(1+EXOG(4,:))./(1-alpha)).^(-epsilon/(1+alpha*(epsilon-1)))); 
    y_vfi(iak_state,:)  = (za*exp(EXOG(2,:)).*EXOG(1,:)).*((k_grid(ik)+EXOG(5,:)).^(alpha)).*(( l_vfi(iak_state,:) - phil).^(1-alpha));
    T_vfi(iak_state,:)  = wage.*EXOG(4,:).*l_vfi(iak_state,:)+wage.*(1+EXOG(4,:))*phil-(EXOG(3,:)+ delta ).*EXOG(5,:);
    pi_vfi(iak_state,:) = D.*(y_vfi(iak_state,:).^((epsilon-1)/epsilon)) - wage.*(1+EXOG(4,:)).* l_vfi(iak_state,:);
    
end

piplusT     = zeros(na,n_choice); 
c           = zeros(na,n_choice);
V_sub       = zeros(na,n_choice);
collateral = chi0 .*AP_mat + chi1.*(exp(KP_mat)-1)  - KP_mat;


while  VFIcontinue==1   
    time1 = tic;
   
    error_vfi_old=error_vfi;    
    
    for i_state = 1:n_state
        prob_mat = repmat(prob(i_state,:),n_choice,1);
        EXP(:,i_state) = sum( prob_mat(:,:).*V_old(:,:) , 2);
    end
    

    for iak_state = 1:n_choice

        ik = rem(iak_state-1,nk)+1; 
        ia = floor((iak_state- 1)/nk)+1;   
        
        piplusT=repmat( pi_vfi(iak_state,:), n_choice,1)+repmat( T_vfi(iak_state,:), n_choice,1);
        
        c = piplusT - (R_mat + delta )*k_grid(ik) + (1+R_mat)*a_grid(ia) - AP_mat - (psi*(KP_mat - k_grid(ik)).*(KP_mat - k_grid(ik)))/(2*k_grid(ik));
        
        CAPeffective = KP_mat + Q_mat;
        
        V_sub = ((1./realpow(abs(c),(gamma-1)*ones(n_choice,n_state)))-1)/(1-gamma) + beta*EXP;

        V_sub(c<=0 | collateral<0 | CAPeffective<=0) = penalty;

        [V_new(iak_state,:), location(iak_state,:)] = max(V_sub);
        [pol_ind_kp(iak_state,:),pol_ind_ap(iak_state,:)] = ind2sub([nk,na],location(iak_state,:));

    end

    
    %%% HOWARD'S IMPROVEMENT ALGORITHM %%%
    Vh_old = V_new;
    iter_howard = 0; 
   
    while iter_howard <= Nh  

        for iak_state = 1: n_choice
            
            ik = rem(iak_state-1,nk)+1; 
            ia = floor((iak_state- 1)/nk)+1;

            c = pi_vfi(iak_state,:) - (EXOG(3,:) + delta )*k_grid(ik) + (1+EXOG(3,:))*a_grid(ia) + T_vfi(iak_state,:)  - a_grid(pol_ind_ap(iak_state,:)) - (psi*(k_grid(pol_ind_kp(iak_state,:)) - k_grid(ik)).*(k_grid(pol_ind_kp(iak_state,:)) - k_grid(ik)))/(2*k_grid(ik));
       
            Vh_new(iak_state,:) = ((1./realpow(abs(c),(gamma-1)*ones(1,nzt*nzp*nr*ntau*nq)))-1)/(1-gamma) + beta .* sum(prob.*Vh_old(location(iak_state,:),:),2)';
        end

        Vh_new = Vh_new.*indicator + Vh_old.*(1-indicator);
                
        Vh_old = Vh_new; 
        iter_howard = iter_howard + 1;  
    end
   
    %%% HOWARD END %%%
    
    if Nh == -1
        error_vfi = max(abs(V_new(:)-V_old(:)));
        V_old = V_new;
    else
        error_vfi = max(abs(Vh_new(:)-V_old(:)));
        V_old = Vh_new;
        
    end
    
    display('----------------------------')        
    display(error_vfi_old)           
    display(error_vfi)            
    
    iter_vfi  = iter_vfi + 1
 
    if error_vfi<error_tolerance
        VFIcontinue=0;
    end      
    
    
    time_elapsed_loop = toc(time1)/3600
       
end

display('----------------------------')
display('----------------------------')
time_VFI=toc/3600
display('----------------------------')
display('----------------------------')

%take results;
V_final             = zeros(na,nk,nzp,nzt,nr,ntau,nq);
a_prime_val         = zeros(na,nk,nzp,nzt,nr,ntau,nq);
k_prime_val         = zeros(na,nk,nzp,nzt,nr,ntau,nq);
c_val               = zeros(na,nk,nzp,nzt,nr,ntau,nq);
binding_collateral  = zeros(na,nk,nzp,nzt,nr,ntau,nq);
sim_pol_ind_ap      = zeros(na,nk,nzp,nzt,nr,ntau,nq);
sim_pol_ind_kp      = zeros(na,nk,nzp,nzt,nr,ntau,nq);
effective_lambda    = zeros(na,nk,nzp,nzt,nr,ntau,nq);

 for i_state = 1:n_state
     
    izp  = EXOG_ind(1,i_state);
    izt  = EXOG_ind(2,i_state);
    ir   = EXOG_ind(3,i_state);
    itau = EXOG_ind(4,i_state);
    iq   = EXOG_ind(5,i_state);   
    for iak_state = 1: n_choice
        
        [ik,ia] = ind2sub([nk,na],iak_state);
                        
        V_final(ia,ik,izp,izt,ir,itau,iq)            = V_new(iak_state, i_state);
        sim_pol_ind_ap(ia,ik,izp,izt,ir,itau,iq)     = pol_ind_ap(iak_state,i_state);
        sim_pol_ind_kp(ia,ik,izp,izt,ir,itau,iq)     = pol_ind_kp(iak_state,i_state);
        a_prime_val(ia,ik,izp,izt,ir,itau,iq)        = a_grid( pol_ind_ap(iak_state,i_state) );
        k_prime_val(ia,ik,izp,izt,ir,itau,iq)        = k_grid( pol_ind_kp(iak_state,i_state) );
        c_val(ia,ik,izp,izt,ir,itau,iq)              = pi_vfi(iak_state,i_state) - (r_grid(ir) + delta )*k_grid(ik) + (1+r_grid(ir))*a_grid(ia) + T_vfi(iak_state,i_state) - a_prime_val(ia,ik,izp,izt,ir,itau,iq) - (psi*(k_prime_val(ia,ik,izp,izt,ir,itau,iq) - k_grid(ik)).^2)./(2*k_grid(ik));
        if pol_ind_kp(iak_state,i_state) < nk 

            effective_lambda(ia,ik,izp,izt,ir,itau,iq)= chi0 + chi1*(exp(k_grid( pol_ind_kp(iak_state,i_state)+1 ))-1)/a_prime_val(ia,ik,izp,izt,ir,itau,iq); 
            if chi0>10^5||chi1>10^5
                effective_lambda(ia,ik,izp,izt,ir,itau,iq)=10^8;
            end                
            if chi0*a_prime_val(ia,ik,izp,izt,ir,itau,iq) + chi1*(exp(k_grid( pol_ind_kp(iak_state,i_state)+1 ))-1)  - k_grid( pol_ind_kp(iak_state,i_state)+1 ) < 0
                binding_collateral(ia,ik,izp,izt,ir,itau,iq) = 1;
            else
                binding_collateral(ia,ik,izp,izt,ir,itau,iq) = 0;
            end             
            
        elseif pol_ind_kp(iak_state,i_state) == nk 
            binding_collateral(ia,ik,izp,izt,ir,itau,iq) = 0;
        end
        
        if V_final(ia,ik,izp,izt,ir,itau,iq) == penalty
            V_final(ia,ik,izp,izt,ir,itau,iq)              = NaN;
            a_prime_val(ia,ik,izp,izt,ir,itau,iq)          = NaN;
            k_prime_val(ia,ik,izp,izt,ir,itau,iq)          = NaN;
            sim_pol_ind_ap(ia,ik,izp,izt,ir,itau,iq)       = NaN;
            sim_pol_ind_kp(ia,ik,izp,izt,ir,itau,iq)       = NaN;
            c_val(ia,ik,izp,izt,ir,itau,iq)                = NaN;
            binding_collateral(ia,ik,izp,izt,ir,itau,iq)   = NaN;
            effective_lambda(ia,ik,izp,izt,ir,itau,iq)     = NaN;            
        end
        
    end
 end



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 %%                      
                %%% 3. Simulation %%%
 
time_sim = tic;

%number of firms and simulation period;
Nfirms      = 50000;
Tperiods    = 1000;
Tshock      = 800;

%generate samples of exogenous variables;
rv_tau = rand(Nfirms,Tperiods);
rv_q = rand(Nfirms,Tperiods);
rv_zt = rand(Nfirms,Tperiods);
       
sim_i_zp_l = 1*ones(1,round(Nfirms*piL));     
sim_i_zp_h = 2*ones(1,round(Nfirms*(1-piL))); 
sim_i_zp = [sim_i_zp_l,sim_i_zp_h];           
        
sim_i_q = zeros(1,Nfirms);        
for i = 1:Nfirms            
    sim_i_q(i) = mod(i,nq) + 1;
end

sim_i_zp_q = [sim_i_zp;sim_i_q];
zp_tau_q_groups = sub2ind([nq,nzp], sim_i_zp_q(2,:), sim_i_zp_q(1,:));

[zp_tau_q_groups,ind_sort_zp_q] = sort(zp_tau_q_groups);

sim_i_zp  = sim_i_zp(ind_sort_zp_q);
sim_i_q   = sim_i_q(ind_sort_zp_q);
sim_i_zp_q = sim_i_zp_q(:,ind_sort_zp_q);

sim_i_q   = repmat(sim_i_q',1,Tperiods);
        
sim_i_tau = zeros(Nfirms,1);
sim_i_tau(:,1) = round(ntau/2); 
for in = 1:Nfirms
    for it = 2:Tperiods
        [~,sim_i_tau(in,it)] = max(rv_tau(in,it) <= cumsum(tauprob(sim_i_tau(in,it-1),:)));
    end    
end

sim_zp_tau_q_distgrid(1) = 0;
for i = 1:nzp*nq
    sim_zp_tau_q_distgrid(i+1) = sim_zp_tau_q_distgrid(i) + sum(zp_tau_q_groups==i);
end
      
sim_i_zt = zeros(Nfirms,Tperiods);
sim_i_zt(:,1) = round(nzt/2);
for in = 1:Nfirms
    for it = 2:Tperiods
        [~,sim_i_zt(in,it)] = max(rv_zt(in,it) <= cumsum(ztprob(sim_i_zt(in,it-1),:)));
    end
end

sim_i_r   = zeros(Tperiods,1);
sim_i_r(1:Tshock-1)         = 6;  %until 1993;
sim_i_r(Tshock:Tshock+17)   =[5,...  %1994
                              5,...  %1995
                              5,...  %1996
                              4,...  %1997
                              3,...  %1998
                              1,...  %1999
                              2,...  %2000
                              2,...  %2001
                              1,...  %2002
                              1,...  %2003
                              1,...  %2004
                              1,...  %2005
                              2,...  %2006
                              3,...  %2007
                              2,...  %2008
                              3,...  %2009
                              2,...  %2010
                              2];  %2011
sim_i_r(Tshock+18:end)=2;               

% pre-allocation;
sim_i_a    = zeros(Nfirms,Tperiods);
sim_i_k    = zeros(Nfirms,Tperiods);
sim_a_val  = zeros(Nfirms,Tperiods);
sim_k_val  = zeros(Nfirms,Tperiods);
sim_c_val  = zeros(Nfirms,Tperiods);
sim_coll_val = zeros(Nfirms,Tperiods);
sim_effective_lambda_val = zeros(Nfirms,Tperiods);
sim_Psi_val = zeros(Nfirms,Tperiods);

a0_init=min(max(a_grid)-0.10,1.00);
k0_init=0.10;
if chi0>10^5||chi1>10^5
    k0_init=0.40;
end
if nq>1
    k0_init=max(-min(q_grid)+0.10,0.10);
end

[value_a0,index_a0]=min(abs(a_grid-a0_init));
[value_k0,index_k0]=min(abs(k_grid-k0_init));

sim_i_a(:,1) = index_a0;
sim_i_k(:,1) = index_k0;
sim_a_val(:,1) = a_grid(sim_i_a(:,1));
sim_k_val(:,1) = k_grid(sim_i_k(:,1));


%generate paths;
it = 2;
index        = sub2ind( size(sim_pol_ind_ap), sim_i_a(:,it-1),sim_i_k(:,it-1),sim_i_zp(:),sim_i_zt(:,it-1),sim_i_r(it-1)*ones(Nfirms,1),sim_i_tau(:,it-1),sim_i_q(:,it-1) );
for it = 2:Tperiods
    sim_i_a(:,it)       = sim_pol_ind_ap( index );
    sim_i_k(:,it)       = sim_pol_ind_kp( index );
    index               = sub2ind( size(c_val), sim_i_a(:,it),sim_i_k(:,it),sim_i_zp(:),sim_i_zt(:,it),sim_i_r(it)*ones(Nfirms,1),sim_i_tau(:,it),sim_i_q(:,it) );
    sim_a_val(:,it)     = a_grid(sim_i_a(:,it));
    sim_k_val(:,it)     = k_grid(sim_i_k(:,it));
    sim_c_val(:,it)     = c_val( index ); 
    sim_coll_val(:,it)  = binding_collateral( index  ); 
    sim_effective_lambda_val(:,it) = effective_lambda(index);
end
sim_zt_val              = zt_grid(sim_i_zt);
sim_zp_val              = repmat(zp_grid(sim_i_zp)',1,Tperiods);
sim_tau_val             = tau_grid(sim_i_tau);
sim_q_val               = q_grid(sim_i_q);


n_chunk = 50;     
n_zp_q_tau_groups = max(zp_tau_q_groups(:));
for in = 1:n_zp_q_tau_groups
    zp_tau_q_frac_group(in) = sum(zp_tau_q_groups==in)/Nfirms;
end
num_firm_in_chunk = round(zp_tau_q_frac_group.*(Nfirms/n_chunk));
min_chunck_size    = sum(num_firm_in_chunk);
for in = 1:n_chunk
    chunk_indicator(:,in) = num_firm_in_chunk'*in;
end
a = chunk_indicator + repmat(sim_zp_tau_q_distgrid(1:n_zp_q_tau_groups)',1,n_chunk);
b = a(:,1) - num_firm_in_chunk';
chunk_indicator_final = [b,a];
if chunk_indicator_final(end,end)>Nfirms
    chunk_indicator_final(:,end) = [];
    n_chunk = n_chunk - 1;
end
iter_all = 1;
for in = 1:n_chunk
    counter = zeros(1,n_zp_q_tau_groups);
    for igroup = 1:n_zp_q_tau_groups
        while counter(igroup) < num_firm_in_chunk(igroup)
            sim_index_order(iter_all) = chunk_indicator_final(igroup,in) + 1  + counter(igroup);
            counter(igroup) = counter(igroup) + 1;
            iter_all = iter_all + 1;
        end   
    end    
end
   
list = sim_index_order;
sim_i_a         = sim_i_a(list,:);     
sim_i_k         = sim_i_k(list,:);     
sim_a_val       = sim_a_val(list,:);
sim_k_val       = sim_k_val(list,:);
sim_c_val       = sim_c_val(list,:);
sim_coll_val    = sim_coll_val(list,:);
sim_effective_lambda_val    = sim_effective_lambda_val(list,:);
sim_i_zt          = sim_i_zt(list,:);
sim_i_zp          = sim_i_zp(list);
sim_i_tau         = sim_i_tau(list,:);
sim_i_q           = sim_i_q(list,:);

sim_r_val         = r_grid(sim_i_r)';
sim_zt_val        = zt_grid(sim_i_zt); 
sim_zp_val        = repmat(zp_grid(sim_i_zp)',1,Tperiods);
sim_tau_val       = tau_grid(sim_i_tau);
sim_q_val         = q_grid(sim_i_q);
sim_productivity_val    = za*exp(sim_zt_val).*sim_zp_val;
sim_l_val               = phil + ((sim_productivity_val).^((epsilon-1)/(1+alpha*(epsilon-1)))).*((sim_k_val+sim_q_val).^(alpha*(epsilon-1)/(1+alpha*(epsilon-1)))).*(D^(epsilon/(1+alpha*(epsilon-1)))).*(mu^(-epsilon/(1+alpha*(epsilon-1)))).*((wage.*(1+sim_tau_val)./(1-alpha)).^(-epsilon/(1+alpha*(epsilon-1))));
sim_y_val               = (sim_productivity_val).*((sim_k_val+sim_q_val).^(alpha)).*((sim_l_val-phil).^(1-alpha));
sim_p_val               = D*sim_y_val.^(-1/epsilon);
sim_rev_val = sim_p_val.*sim_y_val;
sim_mrpk_val            = (alpha/mu)*sim_p_val.*sim_y_val./sim_k_val;
sim_mrpl_val            = ((1-alpha)/mu)*sim_p_val.*sim_y_val./sim_l_val;
sim_b_val               = sim_k_val-sim_a_val;
sim_leverage_val        = sim_b_val./sim_k_val;
sim_Psi_val= [];
for t=1:Tperiods-1
    sim_Psi_val(:,t)        = (psi/2)*((sim_k_val(:,t+1)-sim_k_val(:,t)).^2) ./ sim_k_val(:,t);
end
sim_Psi_val(:,Tperiods)  = sim_Psi_val(:,Tperiods-1);
    
Nfirms = length(sim_index_order);

std_mrpk(:,1)=std(log(sim_mrpk_val(:,:)));
std_mrpl(:,1)=std(log(sim_mrpl_val(:,:)));
logTFP_efficient(:,1) = log( (sum(sim_productivity_val(:,:).^(epsilon-1))).^(1/(epsilon-1)) );
rev(:,1)=sum(sim_rev_val(:,:));
yy(:,1)=sum(sim_y_val(:,:));
YY(:,1)=(sum(sim_y_val(:,:).^((epsilon-1)/epsilon))).^(epsilon/(epsilon-1));
KK(:,1)=sum(sim_k_val(:,:)+sim_q_val(:,:));
LL(:,1)=sum(sim_l_val(:,:))-phil*Nfirms;    
AA(:,1)=sum(sim_a_val(:,:));
CC(:,1)=sum(sim_c_val(:,:));
logTFP_observed(:,1)= log( YY(:,1)./(((KK(:,1)).^(alpha)).*((LL(:,1)).^(1-alpha)))  );
MIS(:,1) = logTFP_observed(:,1)-logTFP_efficient(:,1);
MEAN_coll=mean(sim_coll_val(:,:));
                
display('----------------------------')
display('----------------------------')
time_elapsed_sim = toc(time_sim)/3600
display('----------------------------')
display('----------------------------')


%create dataset;
    
display('----------------------------')  
display('-- creating and saving dataset ---')  
delete dataset.xlsx   
  
dataset=[]; 

if Nfirms>=10000
    Ndataset_start=1;
    Ndataset_end=10000;   
    Ndataset=Ndataset_end-Ndataset_start+1; 
else   
    Ndataset_start=1;  
    Ndataset_end=Nfirms;
    Ndataset=Nfirms;
end

Tdataset_start=Tshock-10;
Tdataset_end=Tshock+26;
Tdataset=Tdataset_end-Tdataset_start+1; 
    
Tshock_dataset=Tshock*ones(Tdataset*Ndataset,1);    
idn=reshape(repmat((1:Ndataset),Tdataset,1),Tdataset*Ndataset,1);
time_dataset=reshape(repmat((Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);
k_dataset=reshape(sim_k_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
a_dataset=reshape(sim_a_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
zt_dataset=reshape(sim_zt_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
zp_dataset=reshape(sim_zp_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
l_dataset=reshape(sim_l_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
y_dataset=reshape(sim_y_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);   
mrpk_dataset=reshape(sim_mrpk_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);  
rev_dataset=reshape(sim_rev_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);  
constrained_dataset=reshape(sim_coll_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);
mrpl_dataset=reshape(sim_mrpl_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1); 
Psi_dataset=reshape(sim_Psi_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);    
c_dataset=reshape(sim_c_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);
effective_lambda_dataset=reshape(sim_effective_lambda_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);   
q_dataset=reshape(sim_q_val(Ndataset_start:Ndataset_end,Tdataset_start:Tdataset_end)',Tdataset*Ndataset,1);   

%note that these aggregates are based on the whole sample of Nfirms;
r_dataset=reshape(repmat(sim_r_val(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);   
logTFPobserved_dataset=reshape(repmat(logTFP_observed(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);   
MIS_dataset=reshape(repmat(MIS(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);  
RR_dataset=reshape(repmat(rev(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);  
yy_dataset=reshape(repmat(yy(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);     
YY_dataset=reshape(repmat(YY(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);     
KK_dataset=reshape(repmat(KK(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);
LL_dataset=reshape(repmat(LL(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);   
AA_dataset=reshape(repmat(AA(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1); 
CC_dataset=reshape(repmat(CC(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);    
COSTRAINED_dataset=reshape(repmat(MEAN_coll(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);    
STD_MRPK_dataset=reshape(repmat(std_mrpk(Tdataset_start:Tdataset_end),1,Ndataset),Tdataset*Ndataset,1);  
   
dataset=[Tshock_dataset,idn,time_dataset,k_dataset,a_dataset,zt_dataset,zp_dataset,l_dataset,y_dataset,mrpk_dataset,rev_dataset,constrained_dataset,r_dataset,logTFPobserved_dataset,MIS_dataset,RR_dataset,yy_dataset,YY_dataset,KK_dataset,LL_dataset,AA_dataset,mrpl_dataset,Psi_dataset,c_dataset,effective_lambda_dataset,CC_dataset,COSTRAINED_dataset,STD_MRPK_dataset,q_dataset];       
xlswrite('dataset.xlsx',dataset);
display('----------------------------')   
display('----------------------------')
time_elapsed_sim = toc(time_sim)/3600
display('----------------------------')
      


