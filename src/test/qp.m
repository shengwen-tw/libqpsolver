pkg load optim

P = [+1 -1;
     -1 +2]
q = [-2;
     -6]

A = []
b = []

lb = [-1;
      -1];
ub = [+1;
      +1];

x = quadprog(P, q, [], [], A, b, lb, ub)
