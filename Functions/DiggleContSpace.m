function [lambda] = DiggleContSpace(spikes,Constants, buff)

% -------------------------------------------------------------
% Function DiggleContSpace
% inputs: coordinates of events, Constants structure
% outputs: Nonparametric estimate of spatial intensity
% Note: This function is overly complicated, it would be better to 
% use standard kernels and then digitise at the end.
%---------------------------------------------------------------

% Initialise
lambda = zeros(Constants.J,Constants.J);
frame = lambda;

% A simple (quick) cylindrical estimator on a digitised map
rabs = buff*Constants.ds;
r = buff;

% Find pixel where event lands 
for i = 1:size(spikes,1)
    %Find x pixel
    [~,x] = find((spikes(i,1) - Constants.s).^2 == min((spikes(i,1) - Constants.s).^2));
    %Find y pixel
    [~,y] = find((spikes(i,2) - Constants.s).^2 == min((spikes(i,2) - Constants.s).^2));
    frame(y,x) = frame(y,x)+1;
end

% Add a digitised (normalised) cylinder for every event in pixel
for i = buff:Constants.J-buff
    for j = buff:Constants.J-buff
        lambda(j,i) = pixincircle(frame,i,j,r)/(pi*rabs^2)/Constants.dt;
    end
end
end