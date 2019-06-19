% Sensitivity analysis to different replicates of birth rate measurements
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
rep2_flow       = rep2_flow(:,:,tpoints2Id); 

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


% Plot BIRTH DATA AND FITS TO ALL DATA AND INDIVIDUAL REPLICATES
figure(1); hold on;
xlabel('Drug concentration (\muM)','FontSize',14,'fontWeight','bold')
ylabel('Birth rates','FontSize',14,'fontWeight','bold')
 
scatter(conc_mean, rep_1_birth,100,'g*', 'LineWidth',2 );
scatter(conc_mean, rep_2_birth,100,'r*', 'LineWidth',2 );
scatter(conc_mean, rep_3_birth,100,'y*', 'LineWidth',2 );

% Generate a range of *birth* rates based on "worst" case scenario where birth rate does not change with drug concentration (tested for reps 1 and 3) to birth rates decreasing with birth concentration according to negative birth rate in paper (from fitting all data)
linetypes   = {'-', '--', 'o', '*',':','-.'}

%% FIT BIRTH DATA
birthData = [rep_1_birth; rep_2_birth; rep_3_birth]

% Fit to each replicate
for i = 1:size(birthData,1)

	pReal       = polyfit(conc_mean' , birthData(i,:), 1);
	if(i==1 || i==3) % Reps 1 and 3 have been tested and the slope is not significantly different from zero. 
		pReal(1)=0
	end
	sens_b  	= polyval(pReal, conc_mean') 
	plot(conc_mean', sens_b, ['k',linetypes{i+1}] ,'LineWidth',2 )
	hold on
	
	save(fullfile(saveDir, ['Birth_rep', num2str(i), '.mat']),'pReal')
end

legend('Replicate 1', 'Replicate 2','Replicate 3', 'Fit Rep 1', 'Fit Rep 2', 'Fit Rep 3','Location', 'southwest')

bfileName   ='Birth_sens'

saveas(gca,fullfile(saveDir,bfileName),'fig')
saveas(gca,fullfile(saveDir,bfileName),'epsc')
saveas(gca,fullfile(saveDir,bfileName),'png')

% PLOT DEATH DATA

deathData       = [rep_1_death; rep_2_death; rep_3_death];

figure(2); hold on;
xlabel('Drug concentration (\muM)','FontSize',14,'fontWeight','bold')
ylabel('death rates','FontSize',14,'fontWeight','bold')
scatter(conc_mean, rep_1_death,100,'g*' , 'LineWidth',2);
scatter(conc_mean, rep_2_death,100 ,'r*', 'LineWidth',2);
scatter(conc_mean, rep_3_death, 100 ,'y*', 'LineWidth',2);
% legend('Replicate 1', 'Replicate 2','Replicate 3','Location', 'northwest')
ylim([0 .035])

% Generate a range of *death* rates 
f3      	= @(b,conc) b(1)*log(b(2)+b(3)*conc);
linecolors 	= ['g', 'r', 'y'];
simRange    = linspace(conc_mean(1), conc_mean(end),100)
eps     	= 1e-6; 

dfileName    ='Death_sens'

% With fit to all data
dData   	 = reshape(deathData', 1, size(deathData,1)*size(deathData,2));
aReal        = lsqcurvefit(f3, [1 1 5], conc', dData, [-Inf eps eps],[Inf Inf Inf]);
modelFit     = f3(aReal,simRange); 
plot(simRange,modelFit, 'k', 'LineWidth',2 )
save(fullfile(saveDir, ['Death_repAll.mat']),'f3','aReal')

sensRange = 1:2
for i = sensRange
	% Simulate a range of fits by increasing by x% the model variable aReal(3) 
    bReal        = [aReal(1), aReal(2), (1+0.5*i)*aReal(3)];
	modelFit     = f3(bReal, simRange); 
	plot(simRange, modelFit, ['k',linetypes{i}])
	hold on
	save(fullfile(saveDir, ['Death_sens', num2str(i), '.mat']),'f3','aReal')
end

sensRange = 1:2
for i = sensRange
	% Simulate a range of fits by increasing by 10% the model variable b(3) 
    bReal        = [aReal(1), aReal(2), (1-0.3*i)*aReal(3)];
	modelFit     = f3(bReal, simRange); 
	plot(simRange, modelFit, ['k',linetypes{i+max(sensRange)}])
	hold on
	save(fullfile(saveDir, ['Death_sens', num2str(i+max(sensRange)), '.mat']),'f3','aReal')
end

legend('Replicate 1', 'Replicate 2','Replicate 3','Fit All', 'Sens 1', 'Sens 2','Sens 3', 'Sens 4','Location', 'northwest')

saveas(gca,fullfile(saveDir,dfileName),'fig')
saveas(gca,fullfile(saveDir,dfileName),'epsc')
saveas(gca,fullfile(saveDir,dfileName),'png')





