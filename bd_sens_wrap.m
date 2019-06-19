% // sensitivity wrapper
% % 
% // using different birth and death rate functions generated from sensitivity_brates
bsens_range = 1:3

dsens_range = 1:6

for j = dsens_range

	for i = bsens_range

		b_sens = ['Birth_rep', num2str(i)]

		d_sens = ['Death_sens', num2str(j)]

		doses_simulation_sens(b_sens, d_sens)

	end
end

