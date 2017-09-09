% Partly based on Giordani & Soderlind (2004). Original code available on https://sites.google.com/site/paulsoderlindecon/home/software

% Modification: uncertainty parameter theta varies over time
% theta_t = 1+e^(a_t)
% a_t = rho*a_(t-1) + (1-rho)*a_bar + std_shock_a*Z,  where Z is standard normal iid

% Clean up
clear all; close all; clc;

%% User input

model = 3;            % 3 for discretion 
if model == 3;        
  theta_lb = 53;      % lower bound for theta (theta_tilde in paper)
  theta_bar = 77;     % unconditional expectation of theta
  std_theta = 70;     % standard deviation of theta
end

theta_const = 0;      % 0: theta varies, 1: theta is constant
if theta_const
  std_theta = 0;
end

% uncertainty variation
rho = 0.5;            % persistence in AR(1) process (cf. Section 4.5)
% properties of {a_t} follow from properties of {theta_t} (cf. Section 4.5)
a_bar = 0.5*log((theta_bar-theta_lb)^4/(std_theta^2+(theta_bar-theta_lb)^2));          
std_shock_a = ((1-rho^2)*log(std_theta^2/(theta_bar-theta_lb)^2+1))^0.5;
% 95% quantile implied by theta_lb, theta_bar, std_theta
% logninv(0.95,a_bar,(std_shock_a^2/(1-rho^2))^0.5) + theta_lb

t = 5;                % number of periods
theta = zeros(t,1);
theta(1,1) = theta_bar;
a = zeros(t,1);
a(1,1) = a_bar;

% parameters of New Keynesian model
gamm = 0.5;           % coefficient of real interest rate in IS curve
alpha = 0.645;        % coefficient of output gap in Philips curve          
rho1 = 0.5;           % persistence in output shocks
rho2 = 0.5;           % persistence in inflation shock
sig1 = 1;             % std. dev. of output shocks
sig2 = 1;             % std. dev. of inflation shocks
% parameters in planner's loss function
bet  = 0.99;          % discount factor
ly   = 0.5;           % coefficient of output gap 
li   = 0.2;           % coefficient of interest rate

IRF = 0;              % 0: shocks occur each period and are drawn from normal distribution
%                       1: IRF to output shock, ie only output shock occurs in initial period and is equal to one std. dev. of its distribution
%                       2: IRF to inflation shock, ie only inflation shock occurs in initial period and is equal to one std. dev. of its distribution
uy = zeros(2,t);      % initialize output shock matrix
upi = zeros(2,t);     % initialize inflation shock matrix
if IRF == 1;
  uy(1,1) = sig1;
elseif IRF == 2;
  upi(2,1) = sig2;      
else
  uy  = randn(2,t);
  uy(2,:) = 0;
  upi = randn(2,t);
  upi(1,:) = 0;
end

DetErrProbIt = 0;     % 1 to calculate detection error probability
NSim = 10000;         % number of simulations for detection error probabilities
TSim = 244;           % sample size for detection error probability

% matrices to save certain values for RE, approximating, and worst-case model for each period
Fu = zeros(3,2);      % saves decision rules
M = zeros(2,2,3);     % saves M for x_(1t+1) = M*x_(1t)+C_1*epsilon_(t+1)
K = zeros(2,2,3);     % saves N for x_(2t) = N*x_(1t)
x1 = zeros(t+1,2,3);  % saves backward-looking variables x1
Sim = zeros(t,3,3);   % saves y,pi,i


%% New Keynesian model of Section 5 in state space form
%  [ 1  0  0   0       [ e_(1t+1)           [ rho1  0     0   0       [e_(1t)      [ 0            [ sig1  0
%    0  1  0   0    *    e_(2t+1)       =      0   rho2   0   0    *   e_(2t)   +    0  * i_t  +     0   sig2 * [ uy_(t+1)
%    0  0  1 gamm        E_(t)y_(t+1)         -1    0     1   0         y_t        gamm              0    0      upi_(t+1)]
%    0  0  0  bet ]      E_(t)y_(t+1) ]        0   -1  -alpha 1 ]      pi_t ]        0 ]             0    0 ]

%       A0          *    [ x1_(t+1)     =              A1          *   [x1_t    +    B * i_t   +       C      * [ uy_(t+1)
%                          x2_(t+1) ]                                   x2_t]                                   upi_(t+1)]

n1 = 2;
n2 = 2;

A0 = eye(n1+n2); 
A0(3,4) = gamm;
A0(4,4) = bet;

A1 = eye(n1+n2);
A1(1,1) = rho1;
A1(2,2) = rho2;
A1(3,1) = -1;
A1(4,2) = -1;
A1(4,3) = -alpha;

B1 = zeros(n1+n2,1);
B1(3) = gamm;

C1 = zeros(n1+n2,n1);
C1(1,1) = sig1;
C1(2,2) = sig2;

A = inv(A0)*A1;
B = inv(A0)*B1;
C = inv(A0)*C1;

Q = zeros(n1+n2,n1+n2);
Q(3,3) = ly;      
Q(4,4) = 1;

R = li;        
U = zeros(n1+n2,1); 

%% Simulation

% optimal rule under RE (obtained by setting the concern for robustness very large, here 1e+5)
if model == 3;               
  [M(:,:,1),K(:,:,1),Ma,Ka,Fu(1,:),Fv,J0,J0a] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,eye(n1),zeros(n2,n1),...
                                       1e-6,0,1,1,1,0,1e+6,C,1e+5);
  C1=C(1:n1,:);                   % cuts off bottom n2 rows ->C1 quadratic (necessary since M is 2x2)
end;

for i = 1:t
   
  % each period the planner solves for the optimally robust rule given his 
  % contemporaneous concern for robustness (cf. Section 4.3)
  if i==1 || theta(i)~=theta(i-1) % rule changes only if theta changes
    if model == 3;    
      [M(:,:,3),K(:,:,3),M(:,:,2),K(:,:,2),Fu(3,:),Fv,J0,J0a] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,eye(n1),zeros(n2,n1),...  %M,Ma: 2x2 matrices
                                       1e-6,0,1,1,1,0,1e+6,C,theta(i));
      Fu(2,:) = Fu(3,:);          % same rule for approximating and worst-case model
    end
  end
  
  % calculation of y,pi,i for RE (j=1), approximating model (j=2) and worst-case model (j=3)
  for j=1:3
    x1(i+1,:,j) = M(:,:,j)*x1(i,:,j)' + C1*(upi(:,i)+uy(:,i));    % for discretion: x1 = x1(t+1) = M*x1(t) + C1*e(t+1)  
    x2 = x1(i+1,:,j)*K(:,:,j)';   % for discretion: x2 = x2(t+1) = K*x1(t)
    if model == 2
      interest = x2(3);
    else
    interest = -x2*Fu(j,:)';      % for discretion: u(t) = -Fu*x1(t)
    end
    Sim(i,1:2,j) = x2(1:2);
    Sim(i,3,j) = interest;
  end 

  % planner's uncertainty in next period
  if i<t                           % not fulfilled in last period
    a(i+1) = (1-rho)*a_bar + rho*a(i) + std_shock_a*randn;
    theta(i+1) = theta_lb+exp(a(i+1));
  end
  
end

%% Evaluation of long-run implications, i.e. for large t

% switches of warning that p-value close to 0 or 0.5
warning('off','stats:jbtest:PTooSmall'); warning('off','stats:jbtest:PTooBig');
% testing worst case y for normality
[jb,pvalue,jbstat] = jbtest(Sim(:,1,3),0.05);
[pvalue,jbstat]
% testing worst case y for autocorrelation
[lb,pvalue,lbstat] = lbqtest(Sim(:,1,3));
[pvalue,lbstat]

% variance
%      RE        approximating   worst-case model
[var(Sim(:,1,1)) var(Sim(:,1,2)) var(Sim(:,1,3));...  % y
 var(Sim(:,2,1)) var(Sim(:,2,2)) var(Sim(:,2,3));...  % pi
 var(Sim(:,3,1)) var(Sim(:,3,2)) var(Sim(:,3,3))]     % i
% 5% confidence interval for std. dev. (interval small for high t)
% [mu,std,muci,stdci] = normfit(Sim(:,3,3))

% autocorrelation
%      RE                                  approximating                              worst-case model
%[corr2(Sim(1:end-1,1,1),Sim(2:end,1,1)) corr2(Sim(1:end-1,1,2),Sim(2:end,1,2)) corr2(Sim(1:end-1,1,3),Sim(2:end,1,3));...  % y
% corr2(Sim(1:end-1,2,1),Sim(2:end,2,1)) corr2(Sim(1:end-1,2,2),Sim(2:end,2,2)) corr2(Sim(1:end-1,2,3),Sim(2:end,2,3));...  % pi
% corr2(Sim(1:end-1,3,1),Sim(2:end,3,1)) corr2(Sim(1:end-1,3,2),Sim(2:end,3,2)) corr2(Sim(1:end-1,3,3),Sim(2:end,3,3))]     % i

% fitting logn, muhat != exp(a_bar+0.5*std_shock_a^2), sigmahat^2 != exp(2*a_bar+std_shock_a^2/(1-rho^2))*(exp(std_shock_a^2/(1-rho^2))-1)
% [par,parci] = lognfit(theta)

%% Graphs
vnames = {'Output gap','Inflation','Interest rate'};
xx = (0:t-1);

% to place legend in an empty corner of the graph
if IRF == 2
  pos = 4;
else 
  pos = 1;
end

figure(1);
  for i=1:3
    subplot(2,2,i);
      plot(xx,Sim(:,i,1),'-k');
      hold on
      plot(xx,Sim(:,i,2),':k');
      plot(xx,Sim(:,i,3),'--k');
      hold off
      if i==1
        legend('RE','Approximating','Worst case',pos);
      end
      title(vnames(i));
  end
  %subplot(2,2,4);
  %  plot(xx,theta(1:end));
  %  title('theta');
         
%% Error detection probability

if DetErrProbIt == 1;
  
  if model == 3;
    %theta = [5.8713, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]; % for rho1 = rho2 = 0 
    theta = [ 46.0365, 47, 48, 50, 55, 60, 70, 80, 90, 100, 150, 200, 250, 400]; % for rho1 = rho2 = 0.5
    %theta = [ 497.369, 500, 510, 550, 600, 700, 800, 900, 1000];  % for rho1 = rho2 = 0.8
    Fu    = NaN;
    Fv0   = NaN;
  end
  constant = 0;

  disp('Calculating error detection probabilities: be patient');
  [rho,JcM] = DetErrProb(A,B,Q,R,U,bet,Fu,n1,n2,cutoff,C,TSim,NSim,theta,123,model,Fv0,1e+6,constant);

  figure(3);
  plot(theta,rho);
  title('Detection error probabilities');
end