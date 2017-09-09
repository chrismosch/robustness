function xM = Var1SimPs(A,epsilon,Tbig,x0)
%Var1SimPs   Calculates impulse response function of a VAR(1) system.
%
%             x(t) = A * x(t-1) +  epsilon(t), where x(t) is nx1
%
%  Usage:     xM = Var1SimPs(A,epsilon,Tbig,x0) or
%                = Var1SimPs(A,epsilon,Tbig)
%
%  Input:     A             nxn VAR(1) matrix, see above
%             epsilon       nx1 or 1xn vector of shocks in inital period or nxTbig
%                           matrix with shocks in all periods
%                           - has to be multiplied by covariance matrix in advance (cf: input in test3 is specified as C1*upi)
%             Tbig          scalar, last period to calculate for
%             x0            nx1 or 1xn vector with starting values, optional
%
%  Output:    xM            Tbig x n matrix, impulse response function
%
%
%
%  Paul.Soderlind@unisg.ch, May 2001
%  URL: https://sites.google.com/site/paulsoderlindecon/home/software
%-----------------------------------------------------------------------

n = size(A,1);

if (nargin == 3);                %if x0 not given as input (cf above), it is generated here 
  x0 = zeros(n,1);
end

x1_t_1 = x0(:);                                %starting vector
xM     = zeros(Tbig,n);                        %to put results in
xM(:)  = NaN;
for t = 1:Tbig;                                %loop over time periods
  x1      = A*x1_t_1 + epsilon(:,t);
  xM(t,:) = x1';                               %storing result for a period t in t-th row of xM
  x1_t_1  = x1;
end