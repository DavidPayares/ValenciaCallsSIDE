function [newsize, flag] = swapdim(size0, dim)
%SWAPDIM   Swapping two adjacent dimensions of an array (DIM and DIM+1).
%   Used only when both A and B are multi-block arrays with 2-D blocks.
%   Example: If the size of A is .......... 5(63)
%            NEWSIZE = SWAPIDS(SIZE0, 2) is 5(36)

    newsize = [size0 1]; % Guarantees that dimension DIM+1 exists.
    newsize = newsize([1:dim-1, dim+1, dim, dim+2:end]);
    flag = true;
end