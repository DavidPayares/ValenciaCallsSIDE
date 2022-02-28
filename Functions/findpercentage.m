function w = findpercentage(x,y,s1,s2,a,b)
% Correction for circle not being completely within square boundary.
w = zeros(size(x,1),1);
Npoints = 50;
for i = 1:size(x,1)
   r = hypot(x(i,1) - y(i,1),x(i,2) - y(i,2)); %find radius
   if (r == 0) || (abs(r - a) > b) %If point is itself or "dr" outside kernel BW set it v. high so 1/w = 0; 
       w(i) = 1e6;
   else 
       d1 = min(abs(s1(1) - x(i,1)),abs(s2(1) - x(i,1)));
       d2 = min(abs(s1(2) - x(i,2)),abs(s2(2) - x(i,2)));
       if r < d1 && r < d2 %circle inside box
           w(i) = 1;
       else %Find it numerically
           THETA=linspace(0,2*pi,Npoints);
           RHO=ones(1,Npoints)*r;
           %Find circle coordinates
           [X,Y] = pol2cart(THETA,RHO);
           X=X+x(i,1);
           Y=Y+x(i,2);
           %Count
           w(i) = sum(((X > s1(1)) & (X < s2(1)))&( (Y > s1(2)) & (Y < s2(2))))/Npoints;
           if (w(i) == 0) 
               w(i) = 1/Npoints; 
           end

       end
       %DEBUG
%        circle(x(1,:),r,1000)
%        rectangle('Position',[0,0,s2(1),s2(2)])
   end
end
end