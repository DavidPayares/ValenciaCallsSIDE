function [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
          squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                          sizeval(idA0,idB0, sizeA0,sizeB0)
%SIZEVAL   Evaluation of array size, block size, and IDs
%    Possible values for IDA and IDB:
%        [DA1 DA2], [DB1 DB2]
%        [DA1 DA2], [DB1]
%        [DA1],     [DB1 DB2]
%        [DA1],     [DB1]
%        [DA1 0],   [0 DB1]
%        [0 DA1],   [DB1 0]
%
%    sizeA/B     Equal to sizeA0/B0 if RESHAPE is not needed in MAIN 1
%    shiftC, delC    Variables controlling MAIN 3.
%    sizeisnew   1x2 logical array; activates reshaping of A and B.
%    idA/B       May change only if squashOK ~= 0
%    squashOK    If only A or B is a multi-block array (M-B) and the other
%                is single-block (1-B), it will be rearranged from N-D to
%                2-D. If both A and B are 1-B or M-B arrays, squashOK = 0.
%                If only A (or B) is a M-B array, squashOK = 1 (or 2).
%    sxtimesOK, timesOK, mtimesOK    Flags controlling MAIN 2 (TRUE/FALSE).
%    sumOK       Dimension along which SUM is performed. If SUM is not
%                needed, sumOK = 0.

% Initializing output arguments

    idA = idA0;
    idB = idB0;
     squashOK = 0;
    sxtimesOK = false;
      timesOK = false;
     mtimesOK = false;
        sumOK = 0;
    shiftC = 0;
    delC = 0;

% Checking for gross input errors

    NidA = numel(idA);
    NidB = numel(idB);
    idA1 = idA(1);
    idB1 = idB(1);
    if  NidA>2 || NidB>2 || NidA==0 || NidB==0 || ...
           ~isreal(idA1) ||    ~isreal(idB1)   || ...
        ~isnumeric(idA1) || ~isnumeric(idB1)   || ...
                 0>idA1  ||          0>idB1    || ... % negative 
         idA1~=fix(idA1) ||  idB1~=fix(idB1)   || ... % non-integer
         ~isfinite(idA1) ||  ~isfinite(idB1) % Inf or NaN               
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['Internal-dimension arguments (e.g., [IDA1 IDA2]) must\n', ...
         'contain only one or two non-negative finite integers']);
    end

% Checking Syntaxes containing zeros (4b/c)

    declared_outer = false;
    idA2 = idA(NidA); % It may be IDA1 = IDA2 (1-D block)
    idB2 = idB(NidB);

    if any(idA==0) || any(idB==0)
        
        % "Inner products": C = MULTIPROD(A, B, [0 DA1], [DB1 0])
        if idA1==0 && idA2>0 && idB1>0 && idB2==0
            idA1 = idA2;
            idB2 = idB1;
        % "Outer products": C = MULTIPROD(A, B, [DA1 0], [0 DB1]) 
        elseif idA1>0 && idA2==0 && idB1==0 && idB2>0
            declared_outer = true;
            idA2 = idA1;
            idB1 = idB2;
        else
            error('MULTIPROD:InvalidDimensionArgument', ...
            ['Misused zeros in the internal-dimension arguments\n', ...
            '(see help heads 4b and 4c)']);
        end
        NidA = 1; 
        NidB = 1;
        idA = idA1;
        idB = idB1;

    elseif (NidA==2 && idA2~=idA1+1) || ...  % Non-adjacent IDs
           (NidB==2 && idB2~=idB1+1)
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['If an array contains 2-D blocks, its two internal dimensions', ... 
        'must be adjacent (e.g. IDA2 == IDA1+1)']);
    end

% ESC - Case for which no reshaping is needed (both A and B are scalars)

    scalarA = isequal(sizeA0, [1 1]);
    scalarB = isequal(sizeB0, [1 1]);
    if scalarA && scalarB
        sizeA = sizeA0;
        sizeB = sizeB0;
        sizeisnew = [false false];
        timesOK = true; return
    end

% Computing and checking adjusted sizes
% The lengths of ADJSIZEA and ADJSIZEB must be >= IDA(END) and IDB(END)

    NsA = idA2 - length(sizeA0); % Number of added trailing singletons
    NsB = idB2 - length(sizeB0);
    adjsizeA = [sizeA0 ones(1,NsA)];
    adjsizeB = [sizeB0 ones(1,NsB)];
    extsizeA = adjsizeA([1:idA1-1, idA2+1:end]); % Size of EDs
    extsizeB = adjsizeB([1:idB1-1, idB2+1:end]);
    p = adjsizeA(idA1);
    q = adjsizeA(idA2);
    r = adjsizeB(idB1);
    s = adjsizeB(idB2);    
    scalarsinA = (p==1 && q==1);
    scalarsinB = (r==1 && s==1);
    singleA = all(extsizeA==1);
    singleB = all(extsizeB==1);
    if q~=r && ~scalarsinA && ~scalarsinB && ~declared_outer
       error('MULTIPROD:InnerDimensionsMismatch', ...
             'Inner matrix dimensions must agree.');
    end

% STEP 1/3 - DIMENSION SHIFTING (FIRST STEP OF AX)
%   Pipeline 1 (using TIMES) never needs left, and may need right shifting.
%   Pipeline 2 (using MTIMES) may need left shifting of A and right of B.

    shiftA = 0;
    shiftB = 0;
    diffBA = idB1 - idA1;    
    if scalarA % Do nothing
    elseif singleA && ~scalarsinB, shiftA = -idA1 + 1; %  Left shifting A
    elseif idB1 > idA1,            shiftA = diffBA;    % Right shifting A        
    end    
    if scalarB % Do nothing
    elseif singleB && ~scalarsinA, shiftB = -idB1 + 1; %  Left shifting B
    elseif idA1 > idB1,            shiftB = -diffBA;   % Right shifting B
    end

% STEP 2/3 - SELECTION OF PROPER ENGINE AND BLOCK SIZE ADJUSTMENTS

    addA  = 0; addB  = 0;
    delA  = 0; delB  = 0;
    swapA = 0; swapB = 0;
    idC1 = max(idA1, idB1);
    idC2 = idC1 + 1;
    checktimes = false;

    if (singleA||singleB) &&~scalarsinA &&~scalarsinB % Engine using MTIMES

        if singleA && singleB 
            mtimesOK = true;
            shiftC=idC1-1; % Right shifting C
            idC1=1; idC2=2;
        elseif singleA
            squashOK = 2;
            idB = [idB1, idB1+1] + shiftB;
        else % singleB
            squashOK = 1;
            idA = [idA1, idA1+1] + shiftA;
        end

        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS
            % OK 
        elseif NidA==2        % 2) 2-D BLOCKS BY 1-D BLOCKS
            addB=idB1+1; delC=idC2;
        elseif NidB==2        % 3) 1-D BLOCKS BY 2-D BLOCKS
            addA=idA1; delC=idC1;
        else                  % 4) 1-D BLOCKS BY 1-D BLOCKS
            if declared_outer
                addA=idA1+1; addB=idB1;
            else
                addA=idA1; addB=idB1+1; delC=idC2;
            end
        end    

    else % Engine using TIMES (also used if SCALARA || SCALARB)
        
        sxtimesOK = true;

        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB || ... % scal-by-mat
                (q==1 && r==1)  % vec-by-vec ("outer")
            elseif p==1 && s==1 % vec-by-vec ("inner")
                swapA=idA1; sumOK=idC1; checktimes=true;
            elseif s==1 % mat-by-vec
                swapB=idB1; sumOK=idC2;
            elseif p==1 % vec-by-mat
                swapA=idA1; sumOK=idC1;
            else % mat-by-mat
                addA=idA2+1; addB=idB1; sumOK=idC2; delC=idC2;
            end

        elseif NidA==2 % 2) 2-D BLOCKS BY 1-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addB=idB1; checktimes=true;
            elseif scalarsinA % scal-by-vec
                delA=idA1;
            elseif scalarsinB % mat-by-scal
                addB=idB1;
            elseif p==1 % vec-by-vec ("inner")
                delA=idA1; sumOK=idC1; checktimes=true;
            else % mat-by-vec
                addB=idB1; sumOK=idC2; delC=idC2;
            end

        elseif NidB==2 % 3) 1-D BLOCKS BY 2-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addA=idA1+1; checktimes=true;
            elseif scalarsinB % vec-by-scal
                delB=idB2;
            elseif scalarsinA % scal-by-mat
                addA=idA1+1;
            elseif s==1 % vec-by-vec ("inner")
                delB=idB2; sumOK=idC1; checktimes=true;
            else % vec-by-mat
                addA=idA1+1; sumOK=idC1; delC=idC1;
            end

        else % 4) 1-D BLOCKS BY 1-D BLOCKS

            if scalarA || scalarB
                timesOK=true;                
            elseif declared_outer % vec-by-vec ("outer")
                addA=idA1+1; addB=idB1;
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB % vec-by-scal
            else % vec-by-vec
                sumOK=idC1; checktimes=true;
            end
        end
    end

% STEP 3/3 - Adjusting the size of A and B. The size of C is adjusted
%            later, because it is not known yet.

    [sizeA, sizeisnew(1)] = adjustsize(sizeA0, shiftA, addA, delA, swapA);
    [sizeB, sizeisnew(2)] = adjustsize(sizeB0, shiftB, addB, delB, swapB);

    if checktimes % Faster than calling BBXFUN
        diff = length(sizeB) - length(sizeA);
        if isequal([sizeA ones(1,diff)], [sizeB ones(1,-diff)])
            timesOK = true;
        end
    end
end