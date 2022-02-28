function Basisnew = FilterBasis(spikes, Basis, shift, scale, rad, bgrate, numweeks, limx, limy, limd)

% For each basis see whether it is significantly outside of Afghanistan. If
% it is remove the basis.
i = 1;
while(i <= Basis.nx)
    if ~inpolygon(Basis.mu1(i)*scale(1)+shift(1),Basis.mu2(i)*scale(2)+shift(2), limx', limy')
        if p_poly_dist(Basis.mu1(i)*scale(1)+shift(1),Basis.mu2(i)*scale(2)+shift(2), limx', limy') > limd
            Basis.nx = Basis.nx-1;
            Basis.mu1(i) = [];
            Basis.mu2(i) = [];
            Basis.sigma2(i) = [];
            Basis.phi(:,:,i) = [];
            i = i-1;
        end
    end
    i = i+1;
end


% For each basis see how many events are within its scope. If there are
% only a few corresponding to a background of exp(-3.5) or less then remove basis.
i=1;
while(i <= Basis.nx)
    distances = hypot(Basis.mu1(i) - spikes(:,1),Basis.mu2(i) - spikes(:,2));
    %     pntsinbasis(i) = sum(distances < 1.2);
    pntsinbasis(i) = sum(distances < rad);
    %     if pntsinbasis(i) < 24
    if pntsinbasis(i) < exp(bgrate)*numweeks*pi*rad^2 %Bg.noise*numofweeks*area
        Basis.nx = Basis.nx-1;
        Basis.mu1(i) = [];
        Basis.mu2(i) = [];
        Basis.sigma2(i) = [];
        Basis.phi(:,:,i) = [];
        i = i-1;
    end
    i = i+1;
end
Basisnew = Basis;
end
