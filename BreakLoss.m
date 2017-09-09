function theta = BreakLoss(theta,model,A,B,Q,R,U,bet,n1,n2,cutoff,C)
  if model == 3;    
      [M,N,Ma,Na,Fu,Fv,J0,J0a,theta] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,eye(n1),zeros(n2,n1),...
                                       1e-6,0,1,1,1,0,1e+4,C,theta);
  end