% Test general sampling algorithm for reversibilities 18/07/18
clc,clearvars,close all

% Number of samples
numSamples = 5e4;
thinning   = 1;

% Example 1: Linear mechanism: i) A <-> P
Omega1  = [1,1,1];
DGr1    = -10;                  % [kJ/mol]
RT      = 8.314*298.15/1e3;     % [kJ/mol]
Rev_1   = generalRevSampling(Omega1,DGr1/RT,numSamples,thinning);

% Example 2: 2 linear promiscuous mechanisms without common steps:
% i) A1 <-> P1 & ii) A2 <-> P2
Omega2  = [1,1,1.05,0,0,0;...
          0,0,0,1,1,1];
DGr2    = [-10;-5];             % [kJ/mol]
Rev_2   = generalRevSampling(Omega2,DGr2/RT,numSamples,thinning);

% Example 3: Linear promiscuous mechanism with 1 common step:
% i) A <-> P1 & ii) A <-> P2
Omega3  = [1,1,1.05,0,0;...
          1,0,0,1,1];
DGr3    = [-10;-8];             % [kJ/mol]
Rev_3    = generalRevSampling(Omega3,DGr3/RT,numSamples,thinning);

% Example 4: Linear promiscuous mechanism with 1 common step and 1 inhibitor:
% i) A <-> P1 & ii) A <-> P2
Omega4  = [1,1,1.05,0,0,0;...
          1,0,0,1,1,0;...
          0,0,0,0,0,1];
DGr4    = [-10;-8;0];           % [kJ/mol]
Rev_4   = generalRevSampling(Omega4,DGr4/RT,numSamples,thinning);

% Compare sampling distributions for example 1 (where we know how the distribution looks like - example 2 is also possible)
% Plot Results Example 1
rev1 = randg(1,size(Omega1,2),numSamples,thinning);
rev1 = rev1./repmat(sum(rev1,1),size(Omega1,2),1);
rev1 = exp(rev1*(DGr1/RT));
figure (1)
xbin = linspace(.05,1,15);
for ix = 1:3
    subplot(2,3,ix)
    f = hist(rev1(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1)
        ylabel('Probability density')
    end
    if (ix == 2)
        title(['Correct distribution: \mu = ',num2str(mean(rev1(ix,:))),', \sigma = ',num2str(std(rev1(ix,:)))])
    else
        title(['\mu = ',num2str(mean(rev1(ix,:))),', \sigma = ',num2str(std(rev1(ix,:)))])
    end
    axis([0,1.05,0,1.2*max(y)])
    subplot(2,3,ix+3)
    f = hist(Rev_1(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1)
        ylabel('Probability density')        
    end
    if (ix == 2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_1(ix,:))),', \sigma = ',num2str(std(Rev_1(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_1(ix,:))),', \sigma = ',num2str(std(Rev_1(ix,:)))])
    end
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end

% Plot Results Example 3
figure (2)
for ix = 1:5
    subplot(2,3,ix)
    f = hist(Rev_3(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
       title(['Distribution new method: \mu = ',num2str(mean(Rev_3(ix,:))),', \sigma = ',num2str(std(Rev_3(ix,:)))])
    else
       title(['\mu = ',num2str(mean(Rev_3(ix,:))),', \sigma = ',num2str(std(Rev_3(ix,:)))]) 
    end    
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end

% Plot Results Example 4
figure (3)
for ix = 1:6
    subplot(2,3,ix)
    f = hist(Rev_4(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_4(ix,:))),', \sigma = ',num2str(std(Rev_4(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_4(ix,:))),', \sigma = ',num2str(std(Rev_4(ix,:)))])
    end
    title(['\mu = ',num2str(mean(Rev_4(ix,:))),', \sigma = ',num2str(std(Rev_4(ix,:)))])
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end
