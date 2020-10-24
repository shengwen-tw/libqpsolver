pkg load optim

P = [+1 -1;
     -1 +2]
q = [-2;
     -6]

A = [ 1 1;
     -1 2;
      2 1]
b = [2;
     2;
     3]

A_eq = []
b_eq = []

lb = [-1;
      -1];
ub = [+3.4;
      +3.3];

x = quadprog(P, q, A, b, A_eq, b_eq, lb, ub)
