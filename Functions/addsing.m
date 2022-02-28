function [newsize, flag] = addsing(size0, dim, ns)
%ADDSING   Adding NS singleton dimensions to the size of an array.
%   Warning: NS is assumed to be a positive integer.
%   Example: If the size of A is ..... SIZE0 = [5 9 3]
%            NEWSIZE = ADDSING(SIZE0, 3, 2) is [5 9 1 1 3]

    if dim > length(size0)
        newsize = size0;
        flag = false;
    else 
        newsize = [size0(1:dim-1), ones(1,ns), size0(dim:end)];
        flag = true;
    end
end