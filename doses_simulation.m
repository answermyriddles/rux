function opt_c = doses_simulation()
close all
clear
M                           = 1e6; % Initial Population size
doses                       = [60 270 360 450 540];
[fitStruct, ~]              = conc_fitting();  % Remember function is not defined at the origin, take care of that
PK_cycle                    = 8; % Instead of assuming one constant drug concentration for 24 hours we will recalculate the current concentration in intervals of size Pk_cycle to estimate concentration and growth rate 
tpop                        = zeros(1,length(doses)); %total for each dose tested
plot_symbols                = {'o-', 'x-', 's-', 'd-','.-'};
mymap                       = colormap(gray(2*length(doses))); % I use a larger color range (2*) the range so that I can choose to exclude black and white later on 
figure; hold on

jj= 1; % simulate 1 week

    for ii = 1:length(doses)

        if jj>1
            M = tpop(ii);
        end
        
        dose            = doses(ii); 
        maxDays         = round(tox2016(dose)); %1-7 days     
        maxHours        = maxDays*24;
        ToffHours       = 7*24-maxHours; 
        % Doses were done as : 0,8,24(next day's 0) but the num cycles refer to when we re-calculate concentration based on dose and this it done according to the PK_cycle interval. 
        % In this case, there are 3 calculations to be done per day at times 0 8 16 (24 is already the next day's 0). That's why we use maxHours/PK_cycle.  +1...
        numCycles       = maxHours/PK_cycle; 
    
        cpop            = zeros(1,round(numCycles));
        conc            = zeros(1,round(numCycles)); 

        % There is no cummulative effect of doses after the first dose so
        % we calculate it separately first
        InitConc        = fitStruct(dose,2); % Concentration of first dose (cmax)
        [InitB, InitD]  = bdrates(InitConc);
        [boff, doff]    = bdrates(0); 

        cpop(1)         = M.*exp((InitB-InitD).*PK_cycle); % growth of the population between first and second dosings
        conc(1)         = InitConc;

        conc(2)         = InitConc + fitStruct(dose,8); % cmax + c8
        [InitB2,InitD2] = bdrates(conc(2)); 
        cpop(2)         = cpop(1).*exp((InitB2-InitD2).*PK_cycle);
       
            % Calculation of growth based on cumulative effect of doses
            % based on an exponentially decaying PK curve of concentration in the blood as a function of initial dose after t time (see conc_fitting() below). 
            
            for i = 3:numCycles  
                
                i8 = 1; i16 = 1; i24 = 1; new_dose = 1;
                
                if (rem(i,3)==0)   
                    new_dose = 0; i24 = 0;                    
                elseif (rem(i,3)==1)
                    i8  = 0;
                elseif (rem(i,3)==2)
                    i16 = 0;
                end

                conc(i)     = InitConc*new_dose + fitStruct(dose,8)*i8 + fitStruct(dose,16)*i16 + fitStruct(dose,24)*i24 ;                 
                [bon, don]  = bdrates(conc(i));  
                cpop(i)     = cpop(i-1)*exp((bon-don)*PK_cycle);

            end % i
       
            tpop (ii)       = cpop(i)*exp((boff-doff)*(ToffHours));   % population after the break (however many days off from treatment on that week)
               
            if i == 21 
                plot(linspace(0+7*(jj-1),maxDays+7*(jj-1), i+1), log2([M cpop]), plot_symbols{ii}, 'Color',mymap(ii+2, :), 'linewidth', 1.5, 'MarkerSize',8 , 'MarkerFaceColor', mymap(ii+2, :)) % M is the initial population which starts the same for all treatment groups            
            else
                plot([linspace(0+7*(jj-1),maxDays+7*(jj-1), i+1) 7*jj], log2([M cpop tpop(ii)]), plot_symbols{ii},'Color', mymap(ii+2, :), 'linewidth', 1.5, 'MarkerSize',8, 'MarkerFaceColor', mymap(ii+2, :))

            end
                
    end % ii
     
xlabel('Days')
ylabel('Simulated JAK2V617F cell population (Log_2)')
hleg1=legend(strcat(cellfun(@num2str,  num2cell(doses), 'UniformOutput', false),' mg/kg'), 'Location', 'southwest');
set(hleg1,'position',[.15 .15 .15 .25])
% set(gca,'yscale','log')
% fileName = 'cellPopulationsinaWeekAllDoses2WeeksLog2';
% fileName = 'cellPopulationsinaWeekAllDosesLog';
fileName = 'cellPopulationsinaWeekFiveDoses1WeekLog2';
saveas(gca,fileName,'fig')
saveas(gca,fileName,'png')
saveas(gca,fileName,'eps')

end

% Fit drug concentration as a function of both intial dose and time
function [fitStruct fitGof]= conc_fitting()
doses       =    [45, 180, 450];
times       = [2, 8, 24];
zdatamean   = [.476, 0.815, 0.01; ...
                56, 1.43, 0.01; ...
                104, 56.4, 0.258]; %uM

dosesx      = repmat(doses', size(zdatamean,2),1);
timesy      = repmat(times, size(zdatamean,1),1);
timesy      = reshape(timesy, size(timesy,1)*size(timesy,2), 1);
zdata       = reshape(zdatamean, size(zdatamean,1)*size(zdatamean,2), 1);

% Add data points for 90 and 270 after 2 hours even though we don't have the other time points for these doses
dosesx      = [dosesx; 90; 270];
timesy      = [timesy; 2; 2];
zdata       = [zdata; 7.35; 74];

conc = @(alpha2, beta2, x, y) ...
        alpha2.*x.*exp(-beta2.*y); 

[dummytest, dummy_gof]= fit([dosesx, timesy], zdata, conc,...
    'StartPoint', [.1, .1], ...
    'Lower', [0, -Inf], ...
    'Robust', 'LAR');

fitStruct   = dummytest;
fitGof      = dummy_gof;

% plot(dummytest, [dosesx, timesy], zdata)
% xlabel('doses')
% ylabel('time')
% zlabel('Plasma concentration')
% saveas(gca,'concentrationSurface','fig')
% saveas(gca,'concentrationSurface','epsc')

end

    