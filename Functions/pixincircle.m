function y = pixincircle(A,xc,yc,radius)

% smooth out (xc,yc) based on the count map A using cylinders of radius "radius" pixels 
% (finds how many circles from events enclose this (xc, yc) )

y=1:size(A,1);
x=1:size(A,2);
[X Y]=meshgrid(x,y);

% This assume the circle falls *entirely* inside the image
R2 = (X-xc).^2+(Y-yc).^2;
[c1,c2] = find(R2 < radius^2);
c(1,:) = c1;
c(2,:) = c2;
c = round(c(:,2:end)); % pixels located in circle
Ac = A(sub2ind(size(A),c(1,:),c(2,:))); % extract value
y = sum(Ac);
end