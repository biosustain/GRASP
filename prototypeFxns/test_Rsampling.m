% Test general sampling algorithm for reversibilities 18/07/18
clc,clearvars,close all

% Number of samples
numSamples = 10000;
thinning   = 10;
burn_in    = 1000;

% Example 1: Linear mechanism: i) A <-> P
Omega1  = [1,1,1];
DGr1    = -10;                  % [kJ/mol]
RT      = 8.314*298.15/1e3;     % [kJ/mol]
Rev_1   = generalRevSampling(Omega1,DGr1/RT,numSamples,thinning, burn_in);

% Example 2: 2 linear promiscuous mechanisms without common steps:
% i) A1 <-> P1 & ii) A2 <-> P2
Omega2  = [1,1,1,0,0,0;...
          0,0,0,1,1,1];
DGr2    = [-10;-5];             % [kJ/mol]
Rev_2   = generalRevSampling(Omega2,DGr2/RT,numSamples,thinning, burn_in);

% Example 3: Linear promiscuous mechanism with 1 common step:
% i) A <-> P1 & ii) A <-> P2
Omega3  = [1,1,1,0;...
          1,1,0,1];
DGr3    = [-10;-8];             % [kJ/mol]
Rev_3    = generalRevSampling(Omega3,DGr3/RT,numSamples,thinning, burn_in);

% Example 4: Linear promiscuous mechanism with 1 common step and 1 inhibitor:
% i) A <-> P1 & ii) A <-> P2
Omega4  = [1,1,1,0,0,0;...
          1,0,0,1,1,0;...
          0,0,0,0,0,1];
DGr4    = [-10;-8;0];           % [kJ/mol]
Rev_4   = generalRevSampling(Omega4,DGr4/RT,numSamples,thinning);


% Example 5: AANAT without inhibition
% i) A + B <-> P1 + Q & A + C <-> P2 + Q
Omega5 = [1, 0, 0, 0, 1, 1, 1, 1;...
          1, 1, 1, 1, 1, 0, 0, 0];
DGr5 = [-30; -10]; % Rev_5 is all NaN with [-30; -2];
Rev_5 = generalRevSampling(Omega5,DGr5/RT,numSamples,thinning, burn_in);


% Example 6: AANAT with inhibition
Omega6 = [1, 0, 0, 0, 0, 1, 1, 1, 1;...
          1, 0, 1, 1, 1, 1, 0, 0, 0];
DGr6 = [-52.7925; -0.5728]; % Rev_5 is Inf and 0 with [-30; -2];
Rev_6 = generalRevSampling(Omega6,DGr6/RT,numSamples,thinning, burn_in);

% Example 7: Radnom mechanism
Omega7 = [0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1;
          0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0;
          1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1;
          1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0];
DGr7 = [-20; -20; -20; -20];
Rev_7 = generalRevSampling(Omega7,DGr7/RT,numSamples,thinning, burn_in);


% Example 8: AANAT with inhibition and an extra track
Omega8 = [1, 0, 0, 0, 0, 1, 1, 1, 1;...
          1, 0, 1, 1, 1, 1, 0, 0, 0;...
          0, 1, 0, 0, 0, 0, 0, 0, 0];
DGr8 = [-52.7925; -0.5728; 0]; % Rev_5 is Inf and 0 with [-30; -2];
Rev_8 = generalRevSampling(Omega8,DGr8/RT,numSamples,thinning, burn_in);

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



% Plot Results Example 2
figure (2)
for ix = 1:6
    subplot(2,3,ix)
    f = hist(Rev_2(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
       title(['Distribution new method: \mu = ',num2str(mean(Rev_2(ix,:))),', \sigma = ',num2str(std(Rev_2(ix,:)))])
    else
       title(['\mu = ',num2str(mean(Rev_2(ix,:))),', \sigma = ',num2str(std(Rev_2(ix,:)))]) 
    end    
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end


% Plot Results Example 3
figure (3)
for ix = 1:4
    subplot(2,2,ix)
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
figure (4)
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


% Plot Results Example 7
figure (7)
for ix = 1:11
    subplot(3,5,ix)
    f = hist(Rev_7(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_7(ix,:))),', \sigma = ',num2str(std(Rev_7(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_7(ix,:))),', \sigma = ',num2str(std(Rev_7(ix,:)))])
    end
    title(['\mu = ',num2str(mean(Rev_7(ix,:))),', \sigma = ',num2str(std(Rev_7(ix,:)))])
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end


% Plot Results Example 8
figure (8)
for ix = 1:9
    subplot(2,5,ix)
    f = hist(Rev_8(ix,:),xbin);
    y = f/trapz(xbin,f);
    bar(xbin,y)
    if (ix == 1) || (ix == 4)
        ylabel('Probability density')
    end
    if (ix==2)
        title(['Distribution new method: \mu = ',num2str(mean(Rev_8(ix,:))),', \sigma = ',num2str(std(Rev_8(ix,:)))])
    else
        title(['\mu = ',num2str(mean(Rev_8(ix,:))),', \sigma = ',num2str(std(Rev_8(ix,:)))])
    end
    title(['\mu = ',num2str(mean(Rev_8(ix,:))),', \sigma = ',num2str(std(Rev_8(ix,:)))])
    axis([0,1.05,0,1.2*max(y)])
    xlabel(['R_',num2str(ix)])
end
