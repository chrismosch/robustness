function [M,K,Ma,Ka,Fu,Fv,J0,J0a,theta] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,Vt1,Kt1,...
                                              ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter,...
                                              C,theta )
%DiscAlgR    Solves for robust optimization problem in discretionary case.
%
%            Same economy and loss function as in ComAlgR, but discretionary
%            equilibrium.
%
%
%
%  Usage:    [M,K,Ma,Ka,Fu,Fv,J0,J0a] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,Vt1,Ct1, ...
%                                                 ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter,...
%                                                 C,theta );
%
%
%  Input:    A        nxn matrix, (n=n1+n2)
%            B        nxk matrix
%            Q        nxn matrix, symmetric
%            R        kxk matrix, symmetric
%            U        nxk matrix
%            bet      scalar, discount factor (eg 0.99)
%            n1       scalar, # of predetermined variables
%            n2       scalar, # of forward looking variables
%            cutoff   scalar, max modulus of stable roots, eg. 1.0001
%            Vt1      n1xn1 matrix. Initial guess for the value function
%                     (a positive definite matrix, such as the identity matrix)
%            Kt1      Initial guess for the matrix K. Ex: an n2xn1 zero matrix
%            ConvCrit scalar, convergence criterion (a small number, such as 1e-6)
%            Vweight  scalar, weight of value function in convergence test
%            Fweight  scalar, weight of policy function in convergence test
%            CritLags scalar, typically 1. If set to x, the convergence criterion is
%                     based on comparison of iteration t with iterations t-1,...,t-x
%            step     scalar in (0,1): factor of updating of V and F as
%                     in Vt1 = step*V + (1-step)*Vt1
%            PrintIt  scalar, 1: printing iteration number 2: printing iteration number
%                     and convergence criteria
%            MaxIter  scalar, maximum number of iterations (eg. 10000). NEW (March 2003)
%            C        nxn1 matrix. CC' is the VCV matrix of the errors;
%                     the last n2 rows are zeros.   cf p11: errors only
%                     part of first two equations
%            theta

%
%  Output:   M        n1xn1 matrix, x1(t+1) = M*x1(t) + C1*e(t+1) in the worst
%                     case solution, C1 is the first n1 rows of C
%            K        n2xn1 matrix, x2 = K*x1
%            Ma       n1xn1 matrix, x1(t+1) = Ma*x1(t) + C1*e(t+1) in the approximating
%                     solution
%            Ka       n2xn1 matrix, x2 = Ka*x1, approximating model
%            Fu       kxn1 matrix, u = -Fu*x1
%            Fv       n1xn1 matrix, v = -Fv*x1
%            J0       scalar, loss function value, worst case
%            J0a      scalar, loss function, approximating model
%
%
%  Note:     The Ka matrix (x2 = Ka*x1) in the approximating model is the same as in the
%            worst case.
%
%  Calls on: DiscAlg
%
%
%  Paolo Giordani and Paul.Soderlind@unisg.ch, Sep 2001
%  URL: https://sites.google.com/site/paulsoderlindecon/home/software
%------------------------------------------------------------------------------

Q = (Q + Q')/2;                %to make symmetric, no impact if Q already symmetric
R = (R + R')/2;

k = size(R,1);
n = n1 + n2;

if any(any(C(n1+1:n,:)));
  warning('Non-zero elements found in the last n2 rows of C')
  return
end

Bstar = [ B, C ];
Ustar = [ U, zeros(n1+n2,n1) ];
Rstar = [ R,            zeros(k,n1); ...
           zeros(n1,k),  -theta*eye(n1) ];

[M,K,V,F] = DiscAlg( A,Bstar,Q,Rstar,Ustar,bet,n1,n2,Vt1,Kt1, ... %theta included via Rstar
                     ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter );

% -----------
% not part of original code                 
if any(eig(V)<=0)|| sum(eig(V))<1  % stops if V not positive-definite (or sum of eigenvalues <1, why?) which implies that theta is below breakdown point
  warning(['theta = ', num2str(theta),' is below breakdown point'])
  Ma=NaN; Ka=NaN; Fu=NaN; Fv=NaN; J0=NaN; J0a=NaN;
  theta = inf;                     % for breakdown-point calculation: theta below breakdown point is not feasible and thus loss for fmincon is set to inf
  return
end
% -----------

A11 = A(1:n1,1:n1);
A12 = A(1:n1,(n1+1):n);
A21 = A((n1+1):n,1:n1);
A22 = A((n1+1):n,(n1+1):n);

B1 = B(1:n1,:);
B2 = B(n1+1:n,:);

C1 = C(1:n1,:);
C2 = C(n1+1:n,:);

Fu = F(1:k,:);
Fv = F(k+1:k+n1,:);

Ma  = A11 + A12*K - B1*Fu;
Ka  = K;                    % same K for approximating and worst case model

Ptilde = [ eye(n1); K; -Fu; -Fv ];             %value function, worst case
W      = Ptilde' * [Q, Ustar; Ustar', Rstar] * Ptilde;
vecV   = (eye(n1*n1) - kron(M',bet*M'))\W(:);
V      = reshape(vecV,n1,n1);                 %solves V = W + bet*M'*V*M

if bet == 1;
  J0 = trace(C1'*V*C1*eye(n1));
else
  J0 = bet/(1-bet) * trace(C1'*V*C1*eye(n1));
end

Ptilde = [ eye(n1); Ka; -Fu; zeros(n1,n1) ];             %value function approximating model
W      = Ptilde' * [Q, Ustar; Ustar', Rstar] * Ptilde;
vecV   = (eye(n1*n1) - kron(Ma',bet*Ma'))\W(:);
V      = reshape(vecV,n1,n1);                 %solves V = W + bet*Ma'*V*Ma
if bet == 1;
  J0a = trace(C1'*V*C1*eye(n1));
else
  J0a = bet/(1-bet) * trace(C1'*V*C1*eye(n1));
end