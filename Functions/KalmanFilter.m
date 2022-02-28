function [Estinfo] = KalmanFilter(y,Initinfo,Constants,FieldMatrices,Basis)


% -------------------------------------------------------------
% Function KalmanFilter
% inputs: y (prior), and ST structures.
% outputs: posterior on initial field
% Description: One-step update of the initial condition
% -------------------------------------------------------------

% Initialize
SigmaR = FieldMatrices.R;
xestpost = zeros(Basis.nx,Constants.N);
sigma2estpost = zeros(Basis.nx,Basis.nx,Constants.N);
xestprior = zeros(Basis.nx,Constants.N);
sigma2estprior = zeros(Basis.nx,Basis.nx,Constants.N);

sigma2estprior(:,:,1) = 50*eye(Basis.nx);
xestprior(:,1) = Initinfo(1).xestprior;
xestpost(:,1) = Initinfo(1).xestpost;

Estinfo = Initinfo;
% Assume point observation
C = eye(Basis.nx);

S = C*sigma2estprior(:,:,1)*C' + SigmaR;
K = sigma2estprior(:,:,1)*C'*inv(S); % Kalman gain
xestpost(:,1) = xestprior(:,1) + K*(y(:,1) - C*xestprior(:,1)); % Mean update
Estinfo(1).xestpost = xestpost(:,1); % Return initial condition
end