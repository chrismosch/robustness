function [M,K,V,F] = DiscAlg( A,B,Q,R,U,bet,n1,n2,Vt1,Kt1,...
                              ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter )
%DiscAlg    Solves the LQ problem under discretion, iterating backwards in time.
%
%
%
%  Usage:     [M,K,V,F] = DiscAlg( A,B,Q,R,U,bet,n1,n2,Vt1,Kt1,...
%                                  ConvCrit,Vweight,Fweight,CritLags,step,PrintIt,MaxIter );
%
%  Input:    A        nxn matrix, (n=n1+n2)
%            B        nxk matrix
%            Q        nxn matrix, symmetric
%            R        kxk matrix, symmetric
%            U        nxk matrix
%            bet      scalar, discount factor (eg 0.99)
%            n1       scalar, # of predetermined variables
%            n2       scalar, # of forward looking variables
%             Vt1        n1xn1 matrix: initial guess of value matrix
%             Kt1        n2xn1 matrix, initial guess of K in x2(t)=K*x1(t)
%             ConvCrit   2x1 convergence criteria for abs(V-Vt1)|abs(F-Ft1)
%             Vweight    scalar or n1xn1 matrix with weights for
%                        difference criterion for V
%             Fweight    scalar or (n1+n2)x1 vector with weights for
%                        difference criterion for F
%             CritLags   #lags of CritVar compared with ConvCrit
%             step       scalar in (0,1): factor of updating of V and F as
%                        in Vt1 = step*V + (1-step)*Vt1
%             PrintIt    1: printing iteration number
%                        2: printing  iteration number and convergence criteria
%             MaxIter    scalar, maximum number of iterations (eg. 10000). NEW (March 2003)
%
%  Output:    M        n1x1 matrix, x1(t+1) = M*x1(t) + e(t+1)
%             K        n2xn1 matrix, x2(t)  = K*x1(t)
%             V        n1xn1 matrix, value function is x1(t)'*V*x1(t)
%             F        kxn1 matrix, decision rule is u(t) = -F*x1(t), where
%                      k is number of elements in u(t)
%
%
%
%  Calls on:  DiscAlg2
%
%
%  Paul Söderlind, Paul.Soderlind@unisg.ch, Aug 2000, Mar 2003
%  URL: https://sites.google.com/site/paulsoderlindecon/home/software
%-----------------------------------------------------------------------

Q = (Q + Q')/2;                %to make symmetric
R = (R + R')/2;

n = n1 + n2;
Ft1 = 1000;

Kdiff = 1000*ones(1+CritLags,2);
iter = 1;
while any( max(Kdiff) > (ConvCrit')) && (iter < MaxIter);   %iterations

  [M,K,F,V] = DiscAlg2( A,B,Q,R,U,bet,n1,n2,Kt1,Vt1 ); %solve period t

  % -----------
  % not part of original code
  if any(eig(V)<0) || sum(eig(V))<1  % stops if V not positive-definite (or sum of eigenvalues <1, why?) which implies that theta is below breakdown point
    M=NaN; K=NaN; F=NaN;
    break
  end
  % -----------
  
  Vdiff = max(  max( Vweight.*abs(V-Vt1) )  );        %changes t+1 -> t
  Fdiff = max(  max( Fweight.*abs(F-Ft1) )  );
  Kdiff = [ Kdiff; ...
            Vdiff, Fdiff ];
  Kdiff = Kdiff(2:size(Kdiff,1),:);                   %latest is last

  Vt1 = step*V + (1-step)*Vt1;                        %"downdating"
  Kt1 = K;
  Ft1 = step*F + (1-step)*Ft1;

  if PrintIt == 1;
    disp(iter);
  elseif PrintIt == 2;
    disp([iter, (max(Kdiff))]);
  end

  iter = iter + 1;

end                                %end iterations

if iter >= MaxIter;
  warning('Maximum number of iterations reached');
end
