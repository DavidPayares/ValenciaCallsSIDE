function [phi] = LocalisedKernelPhi_Cont(s1,s2,mu1,mu2,sigma21,sigma22)

% Evaluate CGRBF with centre (mu1,mu2) with stds (sigma21,sigma22) at
% points (s1,s2)

phi = zeros(length(s1),1);
beta1 = sqrt(pi/sigma21);
beta2 = sqrt(pi/sigma22);
l1 = 2*pi/beta1;
l2 = 2*pi/beta2;
slow1 = (mu1 - l1);  %Find limit kernel
shigh1 = (mu1 + l1);  %Find limit kernel
slow2 = (mu2 - l2);  %Find limit kernel
shigh2 = (mu2 + l2);  %Find limit kernel

Delta1 = beta1*abs(s1 - mu1);
Delta2 = beta2*abs(s2 - mu2);
phi(:,1) = ((2*pi - Delta1).*(1 + cos(Delta1)/2) + 3/2*sin(Delta1))./(3*pi).*((2*pi - Delta2).*(1 + cos(Delta2)/2) + 3/2*sin(Delta2))./(3*pi);
phi = phi.*~(s1 < slow1 | s1 > shigh1 | s2 < slow2 | s2 > shigh2);
end