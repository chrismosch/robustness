
% Clean up
clear all; close all; clc;

%% User input

model = 3;            % 3 for discretion 

theta0 = 10000;       % initial theta

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

cutoff = 1.0;         % scalar, max modulus of stable roots, eg. 1.0001 (only relevant for commitment and simple rule)

%% New Keynesian model of Section 5 in state space form
%  [ 1  0  0   0       [ e_(1t+1)           [ rho1  0     0   0       [e_(1t)      [ 0            [ sig1  0
%    0  1  0   0    *    e_(2t+1)       =      0   rho2   0   0    *   e_(2t)   +    0  * i_t  +     0   sig2 * [ uy_(t+1)
%    0  0  1 gamm        E_(t)y_(t+1)         -1    0     1   0         y_t        gamm              0    0      upi_(t+1)]
%    0  0  0  bet ]      E_(t)y_(t+1) ]        0   -1  -alpha 1 ]      pi_t ]        0 ]             0    0 ]

%       A0          *    [ x1_(t+1)     =              A1          *   [x1_t    +    B * i_t   +       C      * [ uy_(t+1)
%                          x2_(t+1) ]                                   x2_t]                                   upi_(t+1)]

n1 = 2;
n2 = 2;

A0=eye(n1+n2);
A0(3,4)=gamm;
A0(4,4)=bet;

A1=eye(n1+n2);
A1(1,1)=rho1;
A1(2,2)=rho2;
A1(3,1)=-1;
A1(4,2)=-1;
A1(4,3)=-alfa;

B1=zeros(n1+n2,1);
B1(3)=gamm;

C1=zeros(n1+n2,n1);
C1(1,1)=sig1;
C1(2,2)=sig2;

A = inv(A0)*A1;
B = inv(A0)*B1;
C = inv(A0)*C1;

Q=zeros(n1+n2,n1+n2);
Q(3,3)=ly;
Q(4,4)=1;

R = li; 
U = zeros(n1+n2,1); 

%% Breakdown point
warning('off','all');  % Switches off warning 'theta below breakdown' 
%                        in DiscAlgR and ComAlgR off while searching for breakdown point 
[theta,fval,exitflag] = fmincon('BreakLoss',theta0,[],[],[],[],0,[],[],optimset,...
                                model,A,B,Q,R,U,bet,n1,n2,cutoff,C)
warning('on','all');
