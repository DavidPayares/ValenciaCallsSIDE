function c = multiprod(a, b, idA, idB)
%MULTIPROD  Multiplying 1-D or 2-D subarrays contained in two N-D arrays.
%   C = MULTIPROD(A,B) is equivalent  to C = MULTIPROD(A,B,[1 2],[1 2])
%   C = MULTIPROD(A,B,[D1 D2]) is eq. to C = MULTIPROD(A,B,[D1 D2],[D1 D2])
%   C = MULTIPROD(A,B,D1) is equival. to C = MULTIPROD(A,B,D1,D1)
%
%   MULTIPROD performs multiple matrix products, with array expansion (AX)
%   enabled. Its first two arguments A and B are "block arrays" of any
%   size, containing one or more 1-D or 2-D subarrays, called "blocks" (*).
%   For instance, a 563 array may be viewed as an array containing five
%   63 blocks. In this case, its size is denoted by 5(63). The 1 or 2
%   adjacent dimensions along which the blocks are contained are called the
%   "internal dimensions" (IDs) of the array ().
%
%   1) 2-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], [DB1 DB2]) contains the products
%         of the PQ matrices in A by the RS matrices in B. [DA1 DA2] are
%         the IDs of A; [DB1 DB2] are the IDs of B.
%
%   2) 2-D by 1-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], DB1) contains the products of the
%         PQ matrices in A by the R-element vectors in B. The latter are
%         considered to be R1 matrices. [DA1 DA2] are the IDs of A; DB1 is
%         the ID of B.
%
%   3) 1-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, DA1, [DB1 DB2]) contains the products of the 
%         Q-element vectors in A by the RS matrices in B. The vectors in A
%         are considered to be 1Q matrices. DA1 is the ID of A; [DB1 DB2]
%         are the IDs of B.
%
%   4) 1-D BY 1-D BLOCK(S) (*)
%      (a) If either SIZE(A, DA1) == 1 or SIZE(B, DB1) == 1, or both,
%             C = MULTIPROD(A, B, DA1, DB1) returns products of scalars by 
%             vectors, or vectors by scalars or scalars by scalars.
%      (b) If SIZE(A, DA1) == SIZE(B, DB1), 
%             C = MULTIPROD(A, B, [0 DA1], [DB1 0]) or 
%             C = MULTIPROD(A, B, DA1, DB1) virtually turns the vectors
%             contained in A and B into 1P and P1 matrices, respectively,
%             then returns their products, similar to scalar products.
%             Namely, C = DOT2(A, B, DA1, DB1) is equivalent to 
%             C = MULTIPROD(CONJ(A), B, [0 DA1], [DB1 0]).
%      (c) Without limitations on the length of the vectors in A and B,
%             C = MULTIPROD(A, B, [DA1 0], [0 DB1]) turns the vectors
%             contained in A and B into P1 and 1Q matrices, respectively,
%             then returns their products, similar to outer products.
%             Namely, C = OUTER(A, B, DA1, DB1) is equivalent to
%             C = MULTIPROD(CONJ(A), B, [DA1 0], [0 DB1]).
%
%   Common constraints for all syntaxes:
%      The external dimensions of A and B must either be identical or 
%      compatible with AX rules. The internal dimensions of each block
%      array must be adjacent (DA2 == DA1 + 1 and DB2 == DB1 + 1 are
%      required). DA1 and DB1 are allowed to be larger than NDIMS(A) and
%      NDIMS(B). In syntaxes 1, 2, and 3, Q == R is required, unless the
%      blocks in A or B are scalars. 
%
%   Array expansion (AX):
%      AX is a powerful generalization to N-D of the concept of scalar
%      expansion. Indeed, A and B may be scalars, vectors, matrices or
%      multi-dimensional arrays. Scalar expansion is the virtual
%      replication or annihilation of a scalar which allows you to combine
%      it, element by element, with an array X of any size (e.g. X+10,
%      X*10, or []-10). Similarly, in MULTIPROD, the purpose of AX is to
%      automatically match the size of the external dimensions (EDs) of A
%      and B, so that block-by-block products can be performed. ED matching
%      is achieved by means of a dimension shift followed by a singleton
%      expansion:
%      1) Dimension shift (see SHIFTDIM).
%            Whenever DA1 ~= DB1, a shift is applied to impose DA1 == DB1.
%            If DA1 > DB1, B is shifted to the right by DA1 - DB1 steps.
%            If DB1 > DA1, A is shifted to the right by DB1 - DA1 steps.
%      2) Singleton expansion (SX).
%            Whenever an ED of either A or B is singleton and the
%            corresponding ED of the other array is not, the mismatch is
%            fixed by virtually replicating the array (or diminishing it to
%            length 0) along that dimension.
% 
%   MULTIPROD is a generalization for N-D arrays of the matrix
%   multiplication function MTIMES, with AX enabled. Vector inner, outer,
%   and cross products generalized for N-D arrays and with AX enabled are
%   performed by DOT2, OUTER, and CROSS2 (MATLAB Central, file #8782).
%   Elementwise multiplications (see TIMES) and other elementwise binary
%   operations with AX enabled are performed by BAXFUN (MATLAB Central,
%   file #23084). Together, these functions make up the ARRAYLAB toolbox.
%
%   Input and output format:
%      The size of the EDs of C is determined by AX. Block size is
%      determined as follows, for each of the above-listed syntaxes:
%      1) C contains PS matrices along IDs MAX([DA1 DA2], [DB1 DB2]).
%      2) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A         PQ  (2-D)     [DA1 DA2]
%         B         R    (1-D)     DB1
%         C (a)     P    (1-D)     MAX(DA1, DB1)
%         C (b)     PQ  (2-D)     MAX([DA1 DA2], [DB1 DB1+1])
%         ----------------------------------------------------
%         (a) The 1-D blocks in B are not scalars (R > 1).
%         (b) The 1-D blocks in B are scalars (R = 1).
%      3) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A           Q  (1-D)     DA1
%         B         RS  (2-D)     [DB1 DB2]
%         C (a)       S  (1-D)     MAX(DA1, DB1)
%         C (b)     RS  (2-D)     MAX([DA1 DA1+1], [DB1 DB2])
%         ----------------------------------------------------
%         (a) The 1-D blocks in A are not scalars (Q > 1).
%         (b) The 1-D blocks in A are scalars (Q = 1).
%      4)     Array     Block size         ID(s)
%         --------------------------------------------------------------
%         (a) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         MAX(P,Q) (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (b) A         P        (1-D)     DA1
%             B         P        (1-D)     DB1
%             C         1        (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (c) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         PQ      (2-D)     MAX([DA1 DA1+1], [DB1 DB1+1])
%         --------------------------------------------------------------
%
%   Terminological notes:
%   (*) 1-D and 2-D blocks are generically referred to as "vectors" and 
%       "matrices", respectively. However, both may be also called
%       scalars if they have a single element. Moreover, matrices with a
%       single row or column (e.g. 13 or 31) may be also called row
%       vectors or column vectors.
%   () Not to be confused with the "inner dimensions" of the two matrices
%       involved in a product X * Y, defined as the 2nd dimension of X and
%       the 1st of Y (DA2 and DB1 in syntaxes 1, 2, 3).
%
%   Examples:
%    1) If  A is .................... a 5(63)2 array,
%       and B is .................... a 5(34)2 array,
%       C = MULTIPROD(A, B, [2 3]) is a 5(64)2 array.
%
%       A single matrix A pre-multiplies each matrix in B
%       If  A is ........................... a (13)    single matrix,
%       and B is ........................... a 10(34) 3-D array,
%       C = MULTIPROD(A, B, [1 2], [3 4]) is a 10(14) 3-D array.
%
%       Each matrix in A pre-multiplies each matrix in B (all possible
%       combinations)
%       If  A is .................... a (63)5   array,
%       and B is .................... a (34)12 array,
%       C = MULTIPROD(A, B, [1 2]) is a (64)52 array.
%
%   2a) If  A is ........................... a 5(63)2 4-D array,
%       and B is ........................... a 5(3)2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5(6)2   3-D array.
%
%   2b) If  A is ........................... a 5(63)2 4-D array,
%       and B is ........................... a 5(1)2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5(63)2 4-D array.
%
%   4a) If both A and B are .................. 5(6)2   3-D arrays,
%       C = MULTIPROD(A, B, 2) is .......... a 5(1)2   3-D array, while
%   4b) C = MULTIPROD(A, B, [2 0], [0 2]) is a 5(66)2 4-D array
%
%   See also DOT2, OUTER, CROSS2, BAXFUN, MULTITRANSP.

% $ Version: 2.1 $
% CODE      by:            Paolo de Leva
%                          (Univ. of Rome, Foro Italico, IT)    2009 Jan 24
%           optimized by:  Paolo de Leva
%                          Jinhui Bai (Georgetown Univ., D.C.)  2009 Jan 24
% COMMENTS  by:            Paolo de Leva                        2009 Feb 24
% OUTPUT    tested by:     Paolo de Leva                        2009 Feb 24
% -------------------------------------------------------------------------

error( nargchk(2, 4, nargin) ); % Allow 2 to 4 input arguments
switch nargin % Setting IDA and/or IDB
    case 2, idA = [1 2]; idB = [1 2];
    case 3, idB = idA;
end

% ESC 1 - Special simple case (both A and B are 2D), solved using C = A * B

     if ndims(a)==2 && ndims(b)==2 && ...
         isequal(idA,[1 2]) && isequal(idB,[1 2])
         c = a * b; return
     end

% MAIN 0 - Checking and evaluating array size, block size, and IDs

     sizeA0 = size(a);
     sizeB0 = size(b);
     [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
     squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                           sizeval(idA,idB, sizeA0,sizeB0);

% MAIN 1 - Applying dimension shift (first step of AX) and 
%          turning both A and B into arrays of either 1-D or 2-D blocks

     if sizeisnew(1), a = reshape(a, sizeA); end    
     if sizeisnew(2), b = reshape(b, sizeB); end

% MAIN 2 - Performing products with or without SX (second step of AX)

     if squashOK % SQUASH + MTIMES (fastest engine)
         c = squash2D_mtimes(a,b, idA,idB, sizeA,sizeB, squashOK); 
     elseif timesOK % TIMES (preferred w.r. to SX + TIMES)
         if sumOK, c = sum(a .* b, sumOK);
         else      c =     a .* b; end
     elseif sxtimesOK % SX + TIMES
         if sumOK, c = sum(bsxfun(@times, a, b), sumOK);
         else      c =     bsxfun(@times, a, b); end
     elseif mtimesOK % MTIMES (rarely used)
         c = a * b;
     end

% MAIN 3 - Reshaping C (by inserting or removing singleton dimensions)

     [sizeC sizeCisnew] = adjustsize(size(c), shiftC, false, delC, false);
     if sizeCisnew, c = reshape(c, sizeC); end
end