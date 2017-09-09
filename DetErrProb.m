function [rho,JcM] = DetErrProb( A,B,Q,R,U,bet,Fu,n1,n2,cutoff,C, ...
                                 T,NSim,theta,randnstate,ModelType,Fv0,MaxIter,constant )
%ErrDetProb    Simulate error detection probabilites for various theta.
%              Model/equilibrium is chosen by a switch
%
%
%  Usage:    [rho,JcM] = ErrDetProb( A,B,Q,R,U,bet,Fu,n1,n2,cutoff,C, ...
%                                    T,NSim,theta,randnstate,ModelType,Fv0,constant );
%
%
%  Input:    A,B,Q,R,U,bet,Fu,n1,n2,cutoff,C: see ComAlgR, ComAlgSR, DiscAlgR, or SimpAlgR
%            T            scalar, simulated sample length, Ex. 142
%            NSim         scalar, # of simulations, eg 250
%            theta        Kx1 or 1xK vector, values to calculate probs for,
%                         Ex. theta = [ 40; 50; 60];
%            randnstate   scalar, seed of random number generator, Ex. 123
%            ModelType    scalar, switch of model to simulate:
%                            1: ComAlgSR
%                            2: ComAlgR
%                            3: DiscAlgR
%                            4: SimpAlgR  (note: requires optimization routine, see SimpAlgROpt)
%            Fv0         matrix, used only for simple rule, see SimpAlgROpt
%            MaxIter     scalar, maximum number of iterations, used only for discretion, see DiscAlgR
%            constant    scalar, 1 if the first element of x is 1
%
%
%  Output:  rho          Kx1 vector, error detection probabilities
%           JcM          Kx2 matrix, constant in loss function, for
%                        [worst case, approximating model]
%
%
%  Calls on: ComAlgSR, ComAlgR, DiscAlgR, SimpAlgROpt, Var1SimPs
%
%
%  Paolo Giordani and Paul Söderlind, Paul.Soderlind@unisg.ch, Sep 2001
%  URL: https://sites.google.com/site/paulsoderlindecon/home/software
%------------------------------------------------------------------------------


Kbig = length(theta);

rho = zeros(Kbig,1);            %to put results in
rho(:) = NaN;
JcM = zeros(Kbig,2);
JcM(:) = NaN;

for k = 1:Kbig;                      %loop over theta

  if ModelType == 1;                 %commitment, stable
    [Mw,Nw,Ma,Na,Fv,J0,J0a] = ComAlgSR( A,B,Q,R,U,bet,cutoff,C,theta(k) );
  elseif ModelType == 2;             %commitment, unstable
    [Mw,Nw,Ma,Na,Fv,J0,J0a] = ComAlgR(  A,B,Q,R,U,bet,n1,n2,cutoff,C,theta(k) );
  elseif ModelType == 3;             %discretion
    [Mw,Nw,Ma,Na,Fu,Fv,J0,J0a] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,eye(n1),zeros(n2,n1),...
                                         1e-6,0,1,0,1,0,MaxIter,C,theta(k) );
    Fv = [ Fv, zeros(n1,n2) ];           %to make Fv n1xn as required below
  elseif ModelType == 4;
    [Mw,Nw,Ma,Na,J0,J0a,Fv] = SimpAlgROpt( Fv0,A,B,Q,R,U,bet,n1,n2,Fu,cutoff,C,theta(k) );
  end

  if constant == 1;
    Mw2    = Mw(2:end,2:end);
    xmeanw = inv( eye(size(Mw2,1))-Mw2 ) * Mw(2:end,1);
    xmeanw = [1;xmeanw];
    Ma2    = Ma(2:end,2:end);
    xmeana = inv( eye(size(Ma2,1))-Ma2 ) * Ma(2:end,1);
    xmeana = [1;xmeana];
  else
    xmeanw = zeros(size(Mw,1),1);
    xmeana = zeros(size(Ma,1),1);
  end

  randn('state',randnstate);

  rarwN = zeros(NSim,2);
  rarwN(:) = NaN;
  for i = 1:NSim;                                %simulate NSim times
    epsilon = randn(T+1,n1);                     %iid N(0,I)
    if ModelType == 1;
      Xw = Var1SimPs(Mw,C*epsilon,T)         + repmat(xmeanw',T,1);    %simulate VAR(1)
      Xa = Var1SimPs(Ma,C*epsilon,T)         + repmat(xmeana',T,1);
    elseif ModelType == 2;
      x1p2  = Var1SimPs(Mw,(epsilon*C')',T)      + repmat(xmeanw',T,1);
      x2up1 = x1p2*Nw';
      Xw    = [ x1p2(:,1:n1), x2up1(:,1:n2) ];
      x1p2  = Var1SimPs(Ma,(epsilon*C')',T)      + repmat(xmeana',T,1);
      x2up1 = x1p2*Na';
      Xa    = [ x1p2(:,1:n1), x2up1(:,1:n2) ];
    elseif ModelType == 3;
      x1 = Var1SimPs(Mw,(epsilon*C(1:n1,:)')',T) + repmat(xmeanw',T,1);
      Xw = [ x1, x1*Nw' ];
      x1 = Var1SimPs(Ma,(epsilon*C(1:n1,:)')',T) + repmat(xmeana',T,1);
      Xa = [ x1, x1*Nw' ];
    elseif ModelType == 4;
      x1 = Var1SimPs(Mw,(epsilon*C(1:n1,:)')',T) + repmat(xmeanw',T,1);
      Xw = [x1, x1*Nw'];
      x1 = Var1SimPs(Ma,(epsilon*C(1:n1,:)')',T) + repmat(xmeana',T,1);
      Xa = [x1, x1*Nw'];
    end
    vw = -Xw*Fv';
    va = -Xa*Fv';
    ra = sum(sum(va.*va,2))/2 - sum(sum(va.*epsilon(2:T+1,:),2));  %approximating model (cf Sargent equation 9.4.3 p260/247)
    rw = sum(sum(vw.*vw,2))/2 + sum(sum(vw.*epsilon(2:T+1,:),2));  %worst case
    rarwN(i,:) = [ (ra<0), (rw<0) ];
  end

  paw    = mean(rarwN)';                            %frequencies < 0
  rho(k) = mean(paw);                               %(pa + pw)/2
  JcM(k,:) = [J0, J0a];                            %loss fn [worst, approx]

  disp(' ');disp('sigma,   theta,   rho ');         % sigma = = -1/theta
  disp( sprintf('%8.4f ',[-1/theta(k), theta(k), rho(k)]) );

end

