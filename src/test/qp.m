pkg load optim

P = [+1 -1;
     -1 +2]
q = [-2;
     -6]

A = []
b = []

lb = [-20;
      -20];
ub = [+9;
      +20];

x = quadprog(P, q, [], [], A, b, lb, ub)
