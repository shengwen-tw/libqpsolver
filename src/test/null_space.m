A = [+1 +4 +0 +1 -3 +2;
     +2 +8 +1 +1 -4 +6;
     -1 -4 -1 +0 +1 -2;
     +1 +4 +0 +1 -3 +1]

N_A = null(A)

%the n-th column of Q is a basis vector of N(A)
%if the corresponging n-th row of R has all zero entries
[Q, R] = qr(A.')
