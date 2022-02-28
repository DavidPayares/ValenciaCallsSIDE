function [Estinfo] = VBEMSmoother_full_nonlinear(spikes,Constants,Initinfo,FieldMatrices,Basis,Thetainfo)

% -------------------------------------------------------------
% Function VBEMFilter_full_nonlinear
% inputs:
% outputs:
% -------------------------------------------------------------

% Initialize
SigmaW = inv(Thetainfo.Meanprecmat);
beta = Constants.beta;
b = Thetainfo.best;

atm_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.atm,s1,s2);
bank_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.bank ,s1,s2);
bar_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.bar,s1,s2);
cafe_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.cafe,s1,s2);
ind_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.ind,s1,s2);
mark_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.mark ,s1,s2);
night_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.night,s1,s2);
police_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.police ,s1,s2);
pub_cont =  @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.pub,s1,s2);
rest_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.rest ,s1,s2);
taxi_cont =   @(s1,s2) interp2(Constants.s1,Constants.s2,Constants.taxi ,s1,s2);

xestprior = zeros(Basis.nx,Constants.N);
xestpost = zeros(Basis.nx,Constants.N);
xbeta = zeros(Basis.nx,Constants.N);
xestRTS = zeros(Basis.nx,Constants.N);

Sigmaprior = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmapost = zeros(Basis.nx,Basis.nx,Constants.N);
Sigmabeta = zeros(Basis.nx,Basis.nx,Constants.N);
SigmaRTS = zeros(Basis.nx,Basis.nx,Constants.N);

Sigmaprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
Sigmapost(:,:,1) = 30*eye(Basis.nx);
xestpost(:,1) = Initinfo(1).xestpost;

mu = Thetainfo.muest;
theta = Thetainfo.thetaest;
Qinv = inv(SigmaW);
dt = Constants.dt;
Estinfo = Initinfo;

% Options
foptions = [0,1e-4,1e-4,1e-6,0,0,0,0,0,0,0,0,0,0,0,1e-8,0.1,0];
options = foptions;
options(14) = 2000;
% options(9) = 1;
options(2) = 0.1;
options(3) = 0.1;

% Forward pass
for i = 2:Constants.N
    xestprior(:,i) = xestpost(:,i-1) + Constants.dt*theta;
    Sigmaprior(:,:,i) = Sigmapost(:,:,i-1) + SigmaW;
    
    Sigmatilde = inv(inv(Sigmapost(:,:,i-1)) + Qinv);
    Sigmastar = inv(Qinv - Qinv*Sigmatilde*Qinv);
    mustar = Sigmastar*(Qinv*Sigmatilde*(inv(Sigmapost(:,:,i-1))*xestpost(:,i-1) - Qinv*theta*dt) + Qinv*theta*dt);
    Sigmastarinv = inv(Sigmastar);
    
    spikecoords = spikes(i).Coords;
    
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        atm_eval = atm_cont(spikecoords(:,1),spikecoords(:,2));
        bank_eval = bank_cont(spikecoords(:,1),spikecoords(:,2));
        bar_eval = bar_cont(spikecoords(:,1),spikecoords(:,2));
        cafe_eval = cafe_cont(spikecoords(:,1),spikecoords(:,2));
        ind_eval = ind_cont(spikecoords(:,1),spikecoords(:,2));
        mark_eval = mark_cont(spikecoords(:,1),spikecoords(:,2));
        night_eval = night_cont(spikecoords(:,1),spikecoords(:,2));
        police_eval = police_cont(spikecoords(:,1),spikecoords(:,2)); 
        pub_eval = pub_cont(spikecoords(:,1),spikecoords(:,2));
        rest_eval = rest_cont(spikecoords(:,1),spikecoords(:,2));
        taxi_eval = taxi_cont(spikecoords(:,1),spikecoords(:,2)); 
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.sigma2(j),Basis.sigma2(j))';
        end
    else
        phieval = zeros(Basis.nx,size(spikecoords,1));
        atm_eval = 0;
        bank_eval = 0;
        bar_eval = 0;
        cafe_eval = 0;
        ind_eval = 0;
        mark_eval = 0;
        night_eval = 0;
        police_eval = 0; 
        pub_eval = 0;
        rest_eval = 0;
        taxi_eval = 0; 
    end
    
    myint = @(xx)  SingleIntForward(xx,Constants,Basis,Thetainfo);
    myint2 = @(xx)  MultiIntForward(xx,Constants,Basis,Thetainfo);
    f = @(xx) -(sum(mu + b(1)*atm_eval + b(2)*bank_eval + b(3)*bar_eval + b(4)*cafe_eval + b(5)*ind_eval + b(6)*mark_eval+ b(7)*night_eval  + b(8)*police_eval + b(9)*pub_eval + b(10)*rest_eval + b(11)*taxi_eval + beta*phieval'*xx') - myint(xx') - ((xx' - mustar)'/Sigmastar)*(xx' - mustar)./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmastar + mustar'/Sigmastar);
    %[temp,options,t1,t2,t3] = scg(f,xestprior(:,i)',options,gradf);
    [xestpost(:,i)] = scg(f,xestprior(:,i)',options,gradf)';
    myint3 = @(xx) MultiIntForward2(xx,Constants,Basis,Thetainfo);
    temp = inv(Sigmastar) + myint3(xestpost(:,i));
    Sigmapost(:,:,i) = inv(temp);
    [i max(xestpost(:,i)) min(xestpost(:,i))];
    Estinfo(i).xestpost = xestpost(:,i);
    Estinfo(i).xestprior = xestprior(:,i);
    Estinfo(i).PKalman = Sigmapost(:,:,i);
end


% Backward pass
Sigmabeta(:,:,i) = 9*eye(Basis.nx);
xbeta(:,end) = xestpost(:,end);
for i = Constants.N-1:-1:1
    spikecoords = spikes(i+1).Coords;
    if ~isempty(spikecoords)
        phieval = zeros(Basis.nx,size(spikecoords,1));
        atm_eval = atm_cont(spikecoords(:,1),spikecoords(:,2));
        bank_eval = bank_cont(spikecoords(:,1),spikecoords(:,2));
        bar_eval = bar_cont(spikecoords(:,1),spikecoords(:,2));
        cafe_eval = cafe_cont(spikecoords(:,1),spikecoords(:,2));
        ind_eval = ind_cont(spikecoords(:,1),spikecoords(:,2));
        mark_eval = mark_cont(spikecoords(:,1),spikecoords(:,2));
        night_eval = night_cont(spikecoords(:,1),spikecoords(:,2));
        police_eval = police_cont(spikecoords(:,1),spikecoords(:,2)); 
        pub_eval = pub_cont(spikecoords(:,1),spikecoords(:,2));
        rest_eval = rest_cont(spikecoords(:,1),spikecoords(:,2));
        taxi_eval = taxi_cont(spikecoords(:,1),spikecoords(:,2)); 
        for j = 1:Basis.nx
            phieval(j,:) = LocalisedKernelPhi_Cont(spikecoords(:,1),spikecoords(:,2),Basis.mu1(j),Basis.mu2(j),Basis.sigma2(j),Basis.sigma2(j))';
        end
    else
        phieval = zeros(Basis.nx,size(spikecoords,1));
        atm_eval = 0;
        bank_eval = 0;
        bar_eval = 0;
        cafe_eval = 0;
        ind_eval = 0;
        mark_eval = 0;
        night_eval = 0;
        police_eval = 0; 
        pub_eval = 0;
        rest_eval = 0;
        taxi_eval = 0; 
    end
    
    f = @(xx) -(sum(mu + b(1)*atm_eval + b(2)*bank_eval + b(3)*bar_eval + b(4)*cafe_eval + b(5)*ind_eval + b(6)*mark_eval+ b(7)*night_eval  + b(8)*police_eval + b(9)*pub_eval + b(10) * rest_eval + b(11)*taxi_eval + beta*phieval'*xx') - myint(xx') - ((xx' - xbeta(:,i+1))'/Sigmabeta(:,:,i+1))*(xx' - xbeta(:,i+1))./2);
    gradf = @(xx) -(sum(beta*phieval,2)' - myint2(xx')' - xx/Sigmabeta(:,:,i+1) + xbeta(:,i+1)'/Sigmabeta(:,:,i+1));
    mudash = scg(f,xestpost(:,i+1)',options,gradf)';
    Sigmadash = inv(inv(Sigmabeta(:,:,i+1)) + myint3(mudash));
    
    Sigmatilde = inv(inv(Sigmadash) + Qinv);
    Sigmabeta(:,:,i) = inv(Qinv - Qinv*Sigmatilde*Qinv);
    xbeta(:,i) = Sigmabeta(:,:,i)*(-dt*Qinv*theta + Qinv*Sigmatilde*(inv(Sigmadash)*mudash + Qinv*dt*theta));
    
    SigmaRTS(:,:,i) = inv(inv(Sigmapost(:,:,i)) + inv(Sigmabeta(:,:,i)));
    xestRTS(:,i) = SigmaRTS(:,:,i)*(inv(Sigmapost(:,:,i))*xestpost(:,i) + inv(Sigmabeta(:,:,i))*xbeta(:,i));
    
    i
    
    Estinfo(i).xestRTS = xestRTS(:,i);
    Estinfo(i).PRTS = SigmaRTS(:,:,i);
    
    temp = Qinv + inv(Sigmabeta(:,:,i+1))+ myint3(xestRTS(:,i+1));
    Sigmatilde = inv(inv(Sigmapost(:,:,i)) + Qinv);
    Crosscov = Sigmatilde*(Qinv)*inv(temp -  Qinv*Sigmatilde*Qinv);
    
    Estinfo(i).W = xestRTS(:,i)*xestRTS(:,i)' + SigmaRTS(:,:,i);
    Estinfo(i).S = xestRTS(:,i)*xestRTS(:,i+1)' + Crosscov;    
end

%Update posterior for next time step
xestpost(:,1) = xestRTS(:,2) - theta*dt; %initial state to the next iteration
% Sigmapost(:,:,1) = SigmaW*inv(eye(Basis.nx)-A*A'); %initial state
end
