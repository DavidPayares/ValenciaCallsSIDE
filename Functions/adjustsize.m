function [sizeA, sizeisnew] = adjustsize(sizeA0, shiftA, addA, delA, swapA)
% ADJUSTSIZE  Adjusting size of a block array.

    % Dimension shifting (by adding or deleting trailing singleton dim.)
    if     shiftA>0, [sizeA,newA1] = addsing(sizeA0, 1, shiftA);
    elseif shiftA<0, [sizeA,newA1] = delsing(sizeA0, 1,-shiftA); 
    else   sizeA = sizeA0;  newA1  = false;
    end
    % Modifying block size (by adding, deleting, or moving singleton dim.)
    if      addA, [sizeA,newA2] = addsing(sizeA, addA+shiftA, 1); % 1D-->2D 
    elseif  delA, [sizeA,newA2] = delsing(sizeA, delA+shiftA, 1); % 2D-->1D
    elseif swapA, [sizeA,newA2] = swapdim(sizeA,swapA+shiftA); % ID Swapping
    else                 newA2  = false;
    end
    sizeisnew = newA1 || newA2;
end