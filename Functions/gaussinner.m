function y = gaussinner(s1,s2,phi)
% Find the inner product on s1,s2 of a 3-D array of basis functions
nx = size(phi,3);
for i = 1:nx
    for j = 1:nx
        y(i,j) = trapz(s2(:,1),(trapz(s1(1,:),phi(:,:,i).*phi(:,:,j),2)));
    end
end
end