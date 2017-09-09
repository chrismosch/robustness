% Creates histograms of IRFs by estimating IRFs for rolling sample periods of 30 years.

% Clean up
clear all; close all; clc;

%% Set up
% load data
[data, xlstext] = xlsread('data_US.xlsx'); % required format: variables in columns, time in rows
% define label for plots
dates = xlstext(2:end,1);
vnames = xlstext(1,2:end);
% define number of variables and of observations
[nobs, nvar] = size(data);

% parameters set by user
lags = 3;                              % number of lags
leads = 4;                             % number of periods after shock for which IRF is calculated
ma = 30*4;                             % rolling sample specifications: 30y a 4 quarters(cf. De Grauwe Pos1985)
n = nobs-ma+1;                         % remaining sample size i.e. number of rolling means

% matrices to save values
IRF = zeros(nvar,nvar,leads+1,n);      % for IRFs
JBstat = zeros(nvar);                  % for Jarque-Bera test statistic
Pval = zeros(nvar);                    % for p-value of Jerque-Bera test

%% Estimation
% calculate IRF for each rolling sample
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

%% Histograms of IRFs four periods (1 year) after shock

% with summary statistics
figure(1)
for i=1:nvar
    for j=1:nvar
      temp = reshape(IRF(i,j,leads+1,:),[n,1]);        % reshapes IRF into nx1 vector
      % next to each subplot is a pseudo-subplot that contains the summary statistics
      k = (i-1)*nvar*2+(j-1)*2+1;                      % keeps track of current plot position
      subplot(nvar,nvar*2,k)                           % subplots that contain graph
        hist(temp)
        ax = gca; 
        set(ax,'Ygrid','on');
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',[192/255 192/255 192/255]);
        if i==1
          title(strcat('e_{',vnames(j),'}'))
        end
        if j==1
          ylabel(vnames(i))
        end
      subplot(nvar,nvar*2,k+1)                         % subplot for summary statistics
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
    end
end

% only IRF of interest rate to cost-push shock
%{
i = 3;  % index of interest rate
j = 2;  % index of cost-push shock
figure(2)
  temp = reshape(IRF(i,j,leads+1,:),[n,1]);            % reshapes IRF into nx1 vector
  % next to each subplot is a pseudo-subplot that contains the summary statistics
  subplot(1,1*2,1)                                     % subplots that contain graph
    hist(temp)
    ax = gca; 
    set(ax,'Ygrid','on');
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[192/255 192/255 192/255]);
    ylabel('Interest rate')
  subplot(1,1*2,2)                                     % subplot for summary statistics
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
%}

% without summary statistics
%{
figure(3)
for i=1:nvar
    for j=1:nvar
      temp = reshape(IRF(i,j,leads+1,:),[n,1]);        % reshapes IRF into nx1 vector
      subplot(nvar,nvar,(i-1)*nvar+j)   
        hist(temp)
        if i==1
          title(strcat('e_{',vnames(j),'}'))
        end
        if j==1
          ylabel(vnames(i))
        end
    end
end
%}