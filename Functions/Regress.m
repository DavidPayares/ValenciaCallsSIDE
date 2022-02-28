function y = Regress(lambda,Basis,mu,beta)

% Regress an intensity map "lambda" onto the selected basis using mean
% square.

% Initialize
phivec = (reshape(Basis.phi,[],Basis.nx))';
y = zeros(Basis.nx,size(lambda,3));
lambda = lambda + 0.001;  % Some regularisation

% Now find the LS estimate of the initial field
for i = 1:size(lambda,3)
    lambdavec = reshape(lambda(:,:,i),[],1);
    y(:,i) = inv(phivec*phivec')*phivec*((log(lambdavec) - mu)/beta);
end
end