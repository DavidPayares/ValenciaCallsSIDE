function y = MultiIntForward2(x,Constants,Basis,Thetainfo)
% An integral required in optimisation of VB-Laplace
beta = Constants.beta;
mu = Thetainfo.muest;
b = Thetainfo.best;
varb = Thetainfo.varb;
varmu = Thetainfo.varmu;
Uest = Basis.Basisvec*x;
y  = Constants.dt*Constants.ds^2*exp(mu+varmu/2).*(beta^2*Basis.Basisvec'*(Basis.Basisvec.*repmat(exp(b(1)*Constants.atm(:) + varb(1)*Constants.atm(:).^2./2  ...
    + b(2)*Constants.bank(:) + varb(2)*Constants.bank(:).^2./2 ...
    + b(3)*Constants.bar(:) + varb(3)*Constants.bar(:).^2./2 ...
    + b(4)*Constants.cafe(:) + varb(4)*Constants.cafe(:).^2./2 ...
    + b(5)*Constants.ind(:) + varb(5)*Constants.ind(:).^2./2 ...
    + b(6)*Constants.mark(:) + varb(6)*Constants.mark(:).^2./2 ...
    + b(7)*Constants.night(:) + varb(7)*Constants.night(:).^2./2 ...
    + b(8)*Constants.police(:) + varb(8)*Constants.police(:).^2./2 ...
    + b(9)*Constants.pub(:) + varb(9)*Constants.pub(:).^2./2 ...
    + b(10)*Constants.rest(:) + varb(10)*Constants.rest(:).^2./2 ...
    + b(11)*Constants.taxi(:) + varb(11)*Constants.taxi(:).^2./2 ...
    + beta*Uest),1,Basis.nx)));
end