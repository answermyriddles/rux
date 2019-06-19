% called from bd_sens_wrap to run doses_simulation_sens on different birth and death rate functions
function [bon,don] = bdrates_sens(c, b_sens, d_sens)

funcDir = '/Users/delaura/Documents/LevineLab/matlab2016';

    load(fullfile(funcDir, [b_sens, '.mat']))
    load(fullfile(funcDir, [d_sens, '.mat']))

    bon   =   polyval(pReal, c);
	don   =   f3(aReal, c);

end