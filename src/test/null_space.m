A = [16   2   3  13
      5  11  10   8
      9   7   6  12
      4  14  15   1]

null(A)

%the n-th column of Q is a basis vector of N(A)
%if the corresponging n-th row of R has all zero entries
[Q, R] = qr(transpose(A))
