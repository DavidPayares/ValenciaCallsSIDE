function c = squash2D_mtimes(a, b, idA, idB, sizeA, sizeB, squashOK)
% SQUASH2D_MTIMES  Multiproduct with single-block expansion (SBX).
%    Actually, no expansion is performed. The multi-block array is
%    rearranged from N-D to 2-D, then MTIMES is applied, and eventually the
%    result is rearranged back to N-D. No additional memory is required.
%    One and only one of the two arrays must be single-block, and its IDs
%    must be [1 2] (MAIN 1 removes leading singletons). Both arrays
%    must contain 2-D blocks (MAIN 1 expands 1-D blocks to 2-D).

    if squashOK == 1 % A is multi-block, B is single-block (squashing A)

        % STEP 1 - Moving IDA(2) to last dimension
        nd = length(sizeA);
        d2 = idA(2);    
        order = [1:(d2-1) (d2+1):nd d2]; % Partial shifting
        a = permute(a, order); % ...Q

        % STEP 2 - Squashing A from N-D to 2-D  
        q = sizeB(1);
        s = sizeB(2);
        lengthorder = length(order);
        collapsedsize = sizeA(order(1:lengthorder-1)); 
        n = prod(collapsedsize);
        a = reshape(a, [n, q]); % NQ    
        fullsize = [collapsedsize s]; % Size to reshape C back to N-D

    else % B is multi-block, A is single-block (squashing B)

        % STEP 1 - Moving IDB(1) to first dimension
        nd = length(sizeB);
        d1 = idB(1);    
        order = [d1 1:(d1-1) (d1+1):nd]; % Partial shifting
        b = permute(b, order); % Q...

        % STEP 2 - Squashing B from N-D to 2-D  
        p = sizeA(1);
        q = sizeA(2);
        lengthorder = length(order);
        collapsedsize = sizeB(order(2:lengthorder)); 
        n = prod(collapsedsize);
        b = reshape(b, [q, n]); % QN
        fullsize = [p collapsedsize]; % Size to reshape C back to N-D

    end

    % FINAL STEPS - Multiplication, reshape to N-D, inverse permutation
    invorder(order) = 1 : lengthorder;
    c = permute (reshape(a*b, fullsize), invorder);
end