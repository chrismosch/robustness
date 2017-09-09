% Creates a joint graph of estimated and predicted IRFs.
% Roughly, its the code of Sim_IRF

% Clean up
clear all; close all; clc;

%% Empirical IRFs
%% Set up
% load data
[data, xlstext] = xlsread('data_US'); % required format: variables in columns, time in rows
% define label for plots
dates = xlstext(2:end,1);
vnames = xlstext(1,2:end);
% define number of variables and of observations
[nobs, nvar] = size(data);

% parameters set by user
lags = 4;                              % number of lags, determined via Hannan-Quinn IC & Likelihood-ratio test in Emp_VAR_order_selection.m
leads = 4;                             % number of periods after shock for which IRF is calculated
ma = 30*4;                             % rolling mean specifications: 30y a 4 quarters(cf. De Grauwe Pos1985)
n = nobs-ma+1;                         % remaining sample size i.e. number of rolling means

% matrices to save values
IRF = zeros(nvar,nvar,leads+1,n);      % for IRFs
JBstat = zeros(nvar);                  % for Jarque-Bera test statistic
Pval = zeros(nvar);                    % for p-value of Jerque-Bera test

%% Estimation
% calculate IRF for each rolling mean
for i=1:n
    data_i = data(i:i+ma-1,:);
[IRF(:,:,:,i),~,h] = impulseresponse(data_i,1,1:lags,leads,2,0);
end

%% Evaluation
% for ease of comparison with simulated IRF: reorder matrix such that y,pi,i
IRF(:,[2 1],:,:) = IRF(:,[1 2],:,:);
IRF([2 1],:,:,:) = IRF([1 2],:,:,:);
vnames([2 1]) = vnames([1 2]);

% summary statistics
Mean = mean(IRF,4);
Median = median(IRF,4);
Min = min(IRF,[],4);
Max = max(IRF,[],4);
Range = range(IRF,4);
Std = std(IRF,0,4);
Skew = skewness(IRF,0,4);
Kurt = kurtosis(IRF,0,4);

% testing for normality
% switches of warning that p-value close to 0 or 0.5
warning('off','stats:jbtest:PTooSmall'); warning('off','stats:jbtest:PTooBig');
for i=1:nvar
    for j=1:nvar
      temp = reshape(IRF(i,j,leads+1,:),[n,1]);        % reshapes IRF into nx1 vector
      [jb,Pval(i,j),JBstat(i,j)] = jbtest(temp,0.05);
    end
end

%% Predicted IRFs
%% User input
model = 3;                       % 3 for discretion 
if model == 3;        
  theta_lb = 53;                 % lower bound for theta (theta_tilde in paper)
  theta_bar = 77;                % unconditional expectation of theta
  std_theta = 70;                % standard deviation of theta
end

% uncertainty variation
rho = 0.5;                       % persistence in AR(1) process (cf. Section 4.5)
% properties of {a_t} follow from properties of {theta_t} (cf. Section 4.5)
a_bar = 0.5*log((theta_bar-theta_lb)^4/(std_theta^2+(theta_bar-theta_lb)^2));          
std_shock_a = ((1-rho^2)*log(std_theta^2/(theta_bar-theta_lb)^2+1))^0.5;

Trep = 125;                      % number of times the impulse response to a shock is calculated
Tsim = 5;                        % number of periods after shock for which IRF is calculated
theta = zeros(Tsim*Trep,1);      % does not reset but varies each period
theta(1,1) = theta_bar;
a = zeros(Tsim*Trep,1);
a(1,1) = a_bar;

% parameters of New Keynesian model
gamm = 0.5;                      % coefficient of real interest rate in IS curve
alpha = 0.645;                   % coefficient of output gap in Philips curve          
rho1 = 0.5;                      % persistence in output shocks
rho2 = 0.5;                      % persistence in inflation shock
sig1 = 1;                        % std. dev. of output shocks
sig2 = 1;                        % std. dev. of inflation shocks
% parameters in planner's loss function
bet  = 0.99;                     % discount factor
ly   = 0.5;                      % coefficient of output gap 
li   = 0.2;                      % coefficient of interest rate

shock = 2;                         % 1: output shock
%                                    2: supply shock
uy = zeros(2,Tsim);              % initializes output shock matrix
upi = zeros(2,Tsim);             % initializes inflation shock matrix
if IRF == 1;
  uy(1,1) = sig1;
else 
  upi(2,1) = sig2;      
end

% matrices to save certain values
x1a = zeros(Tsim+1,2,Trep);      % saves backward-looking variables for approximating model for each period 
x1w = zeros(Tsim+1,2,Trep);      % saves backward-looking variables for worst-case model for each period
nvar = 3;                        % 3 variables: y,pi,i
SimR_a = zeros(Tsim,nvar,Trep);  % saves y,pi,i for approximating model for each period
SimR_w = zeros(Tsim,nvar,Trep);  % saves y,pi,i for worst-case model for each period

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

%% Robust solution
for j = 1:Trep
for i = 1:Tsim

  T = (j-1)*Tsim + i;                                       % index of overall period
  
  % optimal robust rules
  if model == 3;    
    [M,N,Ma,Na,Fu,Fv,J0,J0a] = DiscAlgR( A,B,Q,R,U,bet,n1,n2,eye(n1),zeros(n2,n1),...
                                       1e-6,0,1,1,1,0,1e+6,C,theta(T));
  end
  
  if model == 3;
    % calculation of y,pi,i for worst-case model
    C1=C(1:n1,:);                                           % cuts off bottom n2 rows ->C1 quadratic (necessary since M is 2x2)
    x1w(i+1,:,j) = M*x1w(i,:,j)' + C1*(upi(:,i)+uy(:,i));   % values of backward-looking variables, i.e. shocks e1,e2    
    x2        = x1w(i+1,:,j)*N';                            % values of forward-looking variables, i.e. y,pi
    interest  = -x2*Fu(end-1:end)';
    SimR_w(i,1:2,j) = x2;
    SimR_w(i,3,j) = interest;
    % calculation of y,pi,i for approximating model
    x1a(i+1,:,j) = Ma*x1a(i,:,j)' + C1*(upi(:,i)+uy(:,i));  % values of backward-looking variables, i.e. shocks e1,e2      
    x2        = x1a(i+1,:,j)*N';                            % values of forward-looking variables, i.e. y,pi
    interest  = -x2*Fu(end-1:end)';
    SimR_a(i,1:2,j) = x2;
    SimR_a(i,3,j) = interest;
  end
  
  % planner's uncertainty in next period
  if i*j<Tsim*Trep                                          % not fulfilled in last period
    a(T+1) = (1-rho)*a_bar + rho*a(T) + std_shock_a*randn;
    theta(T+1) = theta_lb+exp(a(T+1));
  end
  
end
end

%% Evaluation of IRF

% General
% Approximating model
Mean_a = mean(SimR_a,3);
Median_a = median(SimR_a,3);
Min_a = min(SimR_a,[],3);
Max_a = max(SimR_a,[],3);
Range_a = range(SimR_a,3);
Std_a = std(SimR_a,0,3);
Skew_a = skewness(SimR_a,[],3);
Kurt_a = kurtosis(SimR_a,[],3);
% Worst-case model
Mean_w = mean(SimR_w,3);
Median_w = median(SimR_w,3);
Min_w = min(SimR_w,[],3);
Max_w = max(SimR_w,[],3);
Range_w = range(SimR_w,3);
Std_w = std(SimR_w,0,3);
Skew_w = skewness(SimR_w,[],3);
Kurt_w = kurtosis(SimR_w,[],3);

% distribution of IRF values for one period after shock
JBstat_a = zeros(2,3);               % to save Jarque-Bera test statistic for approximating model
Pval_a = zeros(2,3);                 % to save p-value of Jerque-Bera test for approximating model
JBstat_w = zeros(2,3);               %  " for worst-case model
Pval_w = zeros(2,3);                 %  " for worst-case model

% switches of warning that p-value close to 0 or 0.5
warning('off','stats:jbtest:PTooSmall'); warning('off','stats:jbtest:PTooBig');
for i=1:3
  temp_a = reshape(SimR_a(2,i,:),[length(SimR_a),1]);              % reshapes IRF for variable i into (Trep x 1) vector
  [jb,Pval_a(1,i),JBstat_a(1,i)] = jbtest(temp_a,0.05);            % testing for normal distribution
  [jb,Pval_a(2,i),JBstat_a(2,i)] = jbtest(log(abs(temp_a)),0.05);  % testing for lognormal distribution
  temp_w = reshape(SimR_w(2,i,:),[length(SimR_w),1]);
  [jb,Pval_w(1,i),JBstat_w(1,i)] = jbtest(temp_w,0.05);   
  [jb,Pval_w(2,i),JBstat_w(2,i)] = jbtest(log(abs(temp_w)),0.05);
end

%% Graphs
% caption
vnames = {'Output gap','Inflation','Interest rate'};
mtype = {'Empirical','Approximating model','Worst-case model'};

% histograms for IRF for period 1 after shock with summary statistics
j = shock;  % index of empirical shock (1 = demand shock, 2 = cost-push shock, see above)
figure(3);
  for i=1:nvar                                            % for each variable
    % empirical
    subplot(nvar,6,(i-1)*6+1)                             % histogram
      temp = reshape(IRF(i,j,leads+1,:),[n,1]);           % reshapes IRF into nx1 vector
      hist(temp)
      ax = gca; 
      set(ax,'Ygrid','on');
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor',[192/255 192/255 192/255]);
      ylabel(vnames(i))
      if i==1
        title(mtype(1))
      end
    subplot(nvar,6,(i-1)*6+2)                             % summary statistics
    descr = {['Mean:' '          ' num2str(Mean(i,j,leads+1))];...
             ['Median:' '       ' num2str(Median(i,j,leads+1))];...
             ['Min:' '            ' num2str(Min(i,j,leads+1))];...
             ['Max:' '           ' num2str(Max(i,j,leads+1))];...
             ['Range:' '        ' num2str(Range(i,j,leads+1))];...
             ['Std. Dev.:' '    ' num2str(Std(i,j,leads+1))];...
             ['Skewness:' '  ' num2str(Skew(i,j,leads+1))];...
             ['Kurtosis:' '      ' num2str(Kurt(i,j,leads+1))];...
             ['JB:' '              ' num2str(JBstat(i,j))];...
             ['p-value:' '        ' num2str(Pval(i,j))]};
        text(0,.5,descr);
        axis off
    % approximating model
    subplot(nvar,6,(i-1)*6+3)                             % histogram
      temp_a = reshape(SimR_a(2,i,:),[length(SimR_a),1]); % reshapes IRF for variable i into Trep column vector
      hist(temp_a)
      ax = gca; 
      set(ax,'Ygrid','on');
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor',[192/255 192/255 192/255]);
      if i==1
        title(mtype(2))
      end
    subplot(nvar,6,(i-1)*6+4)                             % summary statistics
        descr = {['Mean:' '          ' num2str(Mean_a(2,i))];...
                 ['Median:' '       ' num2str(Median_a(2,i))];...
                 ['Min:' '            ' num2str(Min_a(2,i))];...
                 ['Max:' '           ' num2str(Max_a(2,i))];...
                 ['Range:' '        ' num2str(Range_a(2,i))];...
                 ['Std. Dev.:' '    ' num2str(Std_a(2,i))];...
                 ['Skewness:' '  ' num2str(Skew_a(2,i))];...
                 ['Kurtosis:' '      ' num2str(Kurt_a(2,i))];...
                 ['JB:' '              ' num2str(JBstat_a(1,i))];...
                 ['p-value:' '        ' num2str(Pval_a(1,i))]};
        text(0,.5,descr);
        axis off
    % worst case model
    subplot(nvar,6,(i-1)*6+5)                             % histogram          
      temp_w = reshape(SimR_w(2,i,:),[length(SimR_w),1]); % reshapes IRF for variable i into Trep column vector
      hist(temp_w)
      ax = gca; 
      set(ax,'Ygrid','on');
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor',[192/255 192/255 192/255]);
      if i==1
        title(mtype(3))
      end
    subplot(nvar,6,(i-1)*6+6)                             % summary statistics
        descr = {['Mean:' '          ' num2str(Mean_w(2,i))];...
                 ['Median:' '       ' num2str(Median_w(2,i))];...
                 ['Min:' '            ' num2str(Min_w(2,i))];...
                 ['Max:' '           ' num2str(Max_w(2,i))];...
                 ['Range:' '        ' num2str(Range_w(2,i))];...
                 ['Std. Dev.:' '    ' num2str(Std_w(2,i))];...
                 ['Skewness:' '  ' num2str(Skew_w(2,i))];...
                 ['Kurtosis:' '      ' num2str(Kurt_w(2,i))];...
                 ['JB:' '              ' num2str(JBstat_w(1,i))];...
                 ['p-value:' '        ' num2str(Pval_w(1,i))]};
        text(0,.5,descr);
        axis off
  end