function [newsize, flag] = delsing(size0, dim, ns)
%DELSING   Removing NS singleton dimensions from the size of an array.
%   Warning: Trailing singletons are not removed
%   Example: If the size of A is SIZE0 = [1 1 1 5 9 3]
%            NEWSIZE = DELSING(SIZE, 1, 3) is  [5 9 3]

    if dim > length(size0)-ns % Trailing singletons are not removed
        newsize = size0;
        flag = false;
    else % Trailing singl. added, so NEWSIZE is guaranteed to be 2D or more
        newsize = size0([1:dim-1, dim+ns:end, dim]);
        flag = true;
    end
end