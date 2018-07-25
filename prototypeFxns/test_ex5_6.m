% Test general sampling algorithm for reversibilities 18/07/18
clc,clearvars,close all

% Number of samples
rng('default')          % for reproducibility
numSamples = 5e4;
thinning   = 5e1;
burn_in    = 1e4;
RT      = 8.314*298.15/1e3;     % [kJ/mol]

% Example 5: AANAT without inhibition
% i) A + B <-> P1 + Q & A + C <-> P2 + Q
Omega5 = [1, 0, 0, 0, 1, 1, 1, 1;...
          1, 1, 1, 1, 1, 0, 0, 0];
DGr5  = [-30; -10]; % Rev_5 is all NaN with [-30; -2];
xprev = [];
Rev_5 = generalRevSampling_new(Omega5,DGr5/RT,numSamples,thinning,burn_in,xprev);

% Example 6: AANAT with inhibition
Omega6 = [1, 0, 0, 0, 0, 1, 1, 1, 1;...
          1, 0, 1, 1, 1, 1, 0, 0, 0;...
          0, 1, 0, 0, 0, 0, 0, 0, 0];
DGr6  = [-52.7925; -0.5728; 0];
xprev = [];
Rev_6 = generalRevSampling_new(Omega6,DGr6/RT,numSamples,thinning,burn_in,xprev);

% Compare sampling distributions
nbins = 50;
xbin  = linspace(1/nbins,1,nbins);

% Plot Results Example 5
figure (5)
for ix = 1:8
    subplot(2,4,ix)
    f = hist(Rev_5(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_5(ix,:))),', \sigma = ',num2str(std(Rev_5(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_5(ix,:))),', \sigma = ',num2str(std(Rev_5(ix,:)))])
    end
    title(['\mu = ',num2str(mean(Rev_5(ix,:))),', \sigma = ',num2str(std(Rev_5(ix,:)))])
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end

% Plot Results Example 6
figure (6)
for ix = 1:9
    subplot(3,4,ix)
    f = hist(Rev_6(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_6(ix,:))),', \sigma = ',num2str(std(Rev_6(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_6(ix,:))),', \sigma = ',num2str(std(Rev_6(ix,:)))])
    end
    title(['\mu = ',num2str(mean(Rev_6(ix,:))),', \sigma = ',num2str(std(Rev_6(ix,:)))])
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end