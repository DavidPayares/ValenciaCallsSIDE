function y = kb1(x,b)
% Epanechnikov kernel
y = 3/(4*b)*(1 - x.^2/b^2).*(abs(x(:,1))<b);
end