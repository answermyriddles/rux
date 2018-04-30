% Fitting birth and death rates of cells exposed to varying concentrations of INC424
% with 3 replicate experiments in December 2013 (rep1) and Summer 2014 (rep2&3) 
% KEEP IN MIND that apoptotic population needs to be combined with either
% dead or alive. Flow was done at multiple time points (0, 24 and 48) to 
% capture the correct time range during which cells were affected by drug.
% After the experiments we observed that the correct time range to detect the 
% effect of the drug was 0-48 hrs.
function [pReal, aReal] = fit_nlinfit_birth_death2018_LogDeath() 
close all
clear

saveDir         = '/Users/delaura/Documents/LevineLab/matlab2016'; % legacy

conc_mean=[0 .1 .5 1 2]'; % drug concentrations in micromolar
conc=repmat(conc_mean,3,1); % repeated so we can use for all replicates later

%% REPLICATE 1 (Data from dec 2013 "December flow analysis.pptx")

% % Pick which timepoints to use for growth rate calculation
tpoints1        = [0 48]; 
tpoints1Id      = tpoints1./24 + 1;
% % Flow results (L; D+A)
rep1_flow_0     = [82.6; 4.38+12.2]; % 0hrs  %[L; D+A]
rep1_flow_0     = rep1_flow_0(:,ones(length(conc_mean),1));
% rep1_flow_beg= % 12hrs
rep1_flow_24    = [83.5 81.9 80.2 79.5 78.1; 3.75+12.1 4.01+13.5 3.42+15.8 3.97+15.9 4.15+17]; % 24hrs %[L; D+A]
rep1_flow_48    = [87.4 82.3 73.5 70.9 67.9; 3.43+8.67 3.91+13.1 6.49+19.3 8.24+20.1 9.18+21.8]; % 48hrs  %[L; D+A]
rep1_flow       = cat(3,rep1_flow_0,rep1_flow_24,rep1_flow_48);
rep1_flow       = rep1_flow(:,:,tpoints1Id); 

% % Total cell counts (48hrs was the only one for which we had total # cells/mL, for earlier timepoints I had to calculate based on # viable cells and % viability because K forgot to record the total)
% The viability was taken for 3 different set2 flasks so we have to calc total # cells/ml for each and take the mean
rep1_0hrs_a     = 0.43/.792; % viable: 0.43 mill/mL
rep1_0hrs_b     = 0.72/.818; % viable:, 0.72 mill/mL
rep1_0hrs_c     = 0.75/.75; % viable: 0.75 mill/mL % 0 hrs
rep1_0hrs       = mean([rep1_0hrs_a, rep1_0hrs_b, rep1_0hrs_c]); %  million cells/mL at the begining of experiment

rep1_0hrs       = repmat(rep1_0hrs, 1, length(conc_mean));
rep1_24hrs      = [0.92 0.92 1.12 1.12 1.04]./ ([89.6 91 87.8 87.5 89]./100); % viable cellls/%viable at 24hrs 
rep1_48hrs      = [2.30 2.04 2.12 2.36 2.33]; % we had the total number of cells/mL at 48hrs
rep1            = [rep1_0hrs; rep1_24hrs; rep1_48hrs];

% rep1_vicell_beg_via = repmat(100, 1, length(conc_mean)); % just 100% because I already adjusted the viability above
rep1_vicell     = rep1(tpoints1Id, :);

%% REPLICATE 2 (first replicate from K's experiments in 2014 august; Data from "New Replicate 1 Flow analysis.pptx" )
tpoints2        = tpoints1; 
tpoints2Id      = tpoints1Id;

rep2_flow_0     = [78.1; 2.64+9.2]; % 0hrs
rep2_flow_0     = rep2_flow_0(:,ones(length(conc_mean),1));
rep2_flow_24    = [93.5 93 83.3 77.4 75.4; 2.09+4.16 2.17+4.75 3.54+13 4.24+18.1 4.68+19.7]; %24hrs %[L; D+A]
rep2_flow_48    = [76.5 78.4 50.9 35.3 32.8; 14.3+7.95 9.7+10.4 20.7+25.7 22.8+38.5 30.6+33.8]; %[L; A+D] % 48hrs
rep2_flow       = cat(3,rep2_flow_0,rep2_flow_24,rep2_flow_48);
rep2_flow     = rep2_flow(:,:,tpoints2Id); 

% % Total cell counts (Data from "VICELL READS from both replicates.docx" from Kaitlyn's experiments on August 2014; Kaitlyn only recorded viable cell/ml so we need to calculate the total) 
rep2_0          = 1.39/.93; % total # cells/mL in millions/mL % 0hrs
rep2_0          = rep2_0(:,ones(length(conc_mean),1));
rep2_24         = [1.39 1.82 1.7 1.64 2.04]./([95 93.5 94.9 94 92.6]./100); %24hrs
rep2_48         = [3.11 1.89 1.38 0.9 0.71]./([95.9 91 77.3 67.4 60.4]./100); % 48hrs % cell population
rep2            = [rep2_0; rep2_24; rep2_48];
rep2_vicell     = rep2(tpoints2Id,:);

%% REPLICATE 3 (Second replicate from K's experiments in 2014 august; Data from "New Replicate 2 Flow analysis.pptx" )
tpoints3        = tpoints1; 
tpoints3Id      = tpoints1./24 + 1;

rep3_flow_0     = [92.5; 3.84+3.4]; % 0hrs  %[A+L; D]
rep3_flow_0     = rep3_flow_0(:,ones(length(conc_mean),1));
rep3_flow_24    = [93.8 92.4 86.8 83.3 82.5; 1.86+3.86 1.87+5.26 2.58+10 2.45+13.7 2.88+14.2]; % 24hrs
rep3_flow_48    = [80.2 85.4 76.7 47.4 34.1; 3.46+16.2 6.31+7.79 14.4+7.96 39.5+11.2 51.4+12.6]; % 48hrs  %[A+L; D]
rep3_flow       = cat(3,rep3_flow_0, rep3_flow_24, rep3_flow_48);
rep3_flow       = rep3_flow(:,:,tpoints3Id); 

%% Total cell counts (Data from "VICELL READS from both replicates.docx" from Kaitlyn's experiments on August 2014; Kaitlyn only recorded viable cell/ml so we need to calculate the total) 
rep3_0          = 1.7/.95; % millions/mL
rep3_0          = rep3_0(:,ones(length(conc_mean),1));
rep3_24         = [1.58 1.62 1.63 1.66 1.89]./[.944 .895 .912 .93 .941];
rep3_48         = [3.41 3.23 2.71 2.57 2.35]./([95 91.1 81.8 77.3 73.3]./100); % 48hrs
rep3            = [rep3_0; rep3_24; rep3_48];
rep3_vicell     = rep3(tpoints3Id,:);

% Estimate birth and death rates from flow data for each replicate
[rep_1_birth, rep_1_death]= calculate_birth_death_rates(rep1_flow,rep1_vicell, tpoints1(2)-tpoints1(1));
[rep_2_birth, rep_2_death]= calculate_birth_death_rates(rep2_flow,rep2_vicell, tpoints2(2)-tpoints2(1));
[rep_3_birth, rep_3_death]= calculate_birth_death_rates(rep3_flow,rep3_vicell, tpoints3(2)-tpoints3(1));

%% BIRTH
neps = 1e-2;
rep_1_birth(abs(rep_1_birth)<neps)=0;
rep_2_birth(abs(rep_2_birth)<neps)=0;
rep_3_birth(abs(rep_3_birth)<neps)=0;

% PLOT DATA
figure(3); hold on;
xlabel('Drug concentration (\muM)','FontSize',14,'fontWeight','bold')
ylabel('Birth rates','FontSize',14,'fontWeight','bold')
 
scatter(conc_mean, rep_1_birth,100,'g*', 'LineWidth',2 );
scatter(conc_mean, rep_2_birth,100,'r*', 'LineWidth',2 );
scatter(conc_mean, rep_3_birth,100,'y*', 'LineWidth',2 );
legend('Replicate 1', 'Replicate 2','Replicate 3','Location', 'northwest')

birthData =[rep_1_birth; rep_2_birth; rep_3_birth];

%% FIT BIRTH DATA 
nboot   = 100;
simRange= linspace(0,2.5,100);
bData   = reshape(birthData', 1, size(birthData,1)*size(birthData,2));
simData = zeros(nboot, length(simRange)); 
params  = zeros(nboot,2);

    for ii = 1:nboot
        [bDataSamp, sampIds] = datasample(bData, size(bData,2));
        xSamp                = conc(sampIds);
        P                    = polyfit(xSamp' , bDataSamp, 1);
        params(ii,:)         = P;
        simData(ii,:)        = polyval(P,simRange);
    %     plot(simRange, simData(ii,:)) % plot the result from each fit
    end
 
pReal       = polyfit(conc' , bData, 1);
BirthFit    = polyval(pReal,simRange);
stdSimB     = std(simData,0,1);
s = shadedErrorBar(simRange,BirthFit,stdSimB,[], 2)

set(s.edge,'LineWidth',2)
s.mainLine.LineWidth = 3;

% if we plot the mean overlaid with the fit to the whole data (used as the death function) one can see that they overlap significantly
% plot(simRange,BirthFit, 'k','LineWidth',3 ); hold on
% meanSim  = mean(simData, 1);
% plot(simRange,meanSim,'r-')

set(gca,'FontSize',16,'fontWeight','bold', 'LineWidth', 1)
% ylim([-.05 .05])
ylim([0 .05])

bfileName   ='birthFit_LinFunc';
% saveas(gca,fullfile(saveDir,bfileName),'fig')
% saveas(gca,fullfile(saveDir,bfileName),'epsc')

%% DEATH RATE
% DATA
rep_1_death(abs(rep_1_death)<neps)=0;
rep_2_death(abs(rep_2_death)<neps)=0;
rep_3_death(abs(rep_3_death)<neps)=0;

deathData       = [rep_1_death; rep_2_death; rep_3_death];

% PLOT DATA
figure(4); hold on;
xlabel('Drug concentration (\muM)','FontSize',14,'fontWeight','bold')
ylabel('death rates','FontSize',14,'fontWeight','bold')
scatter(conc_mean, rep_1_death,100,'g*' , 'LineWidth',2);
scatter(conc_mean, rep_2_death,100 ,'r*', 'LineWidth',2);
scatter(conc_mean, rep_3_death, 100 ,'y*', 'LineWidth',2);
legend('Replicate 1', 'Replicate 2','Replicate 3','Location', 'northwest')

% FIT LOG FUNCTION TO DATA
% Bootstrap data and fit multiple times
dData   = reshape(deathData', 1, size(deathData,1)*size(deathData,2));
simData = zeros(nboot, length(simRange)); 
params  = zeros(nboot, 3); 
f3      = @(b,conc) b(1)*log(b(2)+b(3)*conc);
eps     = 1e-6; 

for ii = 1:nboot
    % bootstrap the data
    [dDataSample, idSample] = datasample(dData, size(dData,2)); % sample with replacement
    concSample  = conc(idSample); % take the equivalent sample of the matching concentrations
    a3              = lsqcurvefit(f3, [1 1 5], concSample', dDataSample, [-Inf eps eps],[Inf Inf Inf]);
    params(ii,:)    = a3;
    simData(ii,:)   = f3(a3,simRange);
%     plot(simRange, simData(ii,:)) % plot the result from each fit
end

%% PLOT SIMULATIONS
figure(4) 
aReal        = lsqcurvefit(f3, [1 1 5], conc', dData, [-Inf 1 -Inf],[Inf Inf Inf]);
modelFit     = f3(aReal,simRange); 
stdSimD      = std(simData,0,1); % bootstrapping error
s            = shadedErrorBar(simRange,modelFit,stdSimD,[], 2);

set(s.edge,'LineWidth',2)
s.mainLine.LineWidth = 3;

set(gca,'FontSize',16,'fontWeight','bold', 'LineWidth', 1)
ylim([0 .05])

dfileName   ='deathFit_logFunc';
% saveas(gca,fullfile(saveDir,dfileName),'fig')
% saveas(gca,fullfile(saveDir,dfileName),'epsc')

save(fullfile(saveDir, 'Birth_dlog_w5.mat'),'pReal')
save(fullfile(saveDir, 'Death_log_w5.mat'),'f3','aReal')

%% Plot G = B-D rates

figure(5); hold on 
% Plot data
xlabel('Drug concentration (\muM)','FontSize',14,'fontWeight','bold')
ylabel('Growth rates (cell/hour)','FontSize',14,'fontWeight','bold')
 
scatter(conc_mean, rep_1_birth-rep_1_death,100,'g*', 'LineWidth',2 );
scatter(conc_mean, rep_2_birth-rep_2_death,100,'r*', 'LineWidth',2 );
scatter(conc_mean, rep_3_birth-rep_3_death,100,'y*', 'LineWidth',2 );
legend('Replicate 1', 'Replicate 2','Replicate 3','Location', 'northwest')

% Plot simulation
stdSim      = std(simData,0,1); % bootstrapping error
plot(simRange,BirthFit-modelFit)
s_g = shadedErrorBar(simRange,BirthFit-modelFit,stdSimB-stdSimD,[],2)

set(s.edge,'LineWidth',2)
s_g.mainLine.LineWidth = 3;

set(gca,'FontSize',20,'fontWeight','bold', 'LineWidth', 1)

ylim([-0.005 .04])

dfileName   ='GrowthRate_BminusD_polyDeath2018';
saveas(gca,fullfile(saveDir,dfileName),'fig')
saveas(gca,fullfile(saveDir,dfileName),'epsc')

end

function [birth, death]=calculate_birth_death_rates(flow, vicell, dt)
delta_time = dt;

% vector of total cells (million/mL) for a particular day for every drug concentration. 
total_beg  = vicell(1,:);
total_end  = vicell(2,:);

flow_beg   = flow(:,:,1)./100; % should be a matrix with columns [live apopt dead] 
flow_end   = flow(:,:,2)./100; % percentages and different concentrations on rows

live_beg    = total_beg.*flow_beg(1,:);
live_end    = total_end.*flow_end(1,:);

dead_beg    = total_beg.*flow_beg(2,:); %(millions)
dead_end   = total_end.*flow_end(2,:); %(millions)

% bd is the birth minus death rate
bd          = log(live_end./live_beg)./delta_time;
death_rate  = (dead_end-dead_beg).*bd./(live_end-live_beg);
birth_rate  = bd+death_rate;
birth       = birth_rate;
death       = death_rate;

end

