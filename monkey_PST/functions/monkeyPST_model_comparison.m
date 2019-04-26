function [R2, F, index_best_model, winning, F_day] = monkeyPST_model_comparison(m, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, all_model_architecture)

for mod = 1: length(all_model_architecture)
    model_architecture   = all_model_architecture{mod};
    [R2(mod), F(mod), F_day(:,mod)] = get_R2_F(m, el, d, electrodes, phase, spatial_ref, binned_trials, n_bins, clip, model_architecture);
    disp(['Model ', all_model_architecture{mod}, ': ', 'R^2 = ',num2str(R2(mod)),' & F = ',num2str(F(mod))])
end

index_best_model = find(F==max(F));
index_wo_best_model = find(1:length(all_model_architecture) ~= index_best_model);
index_second_model = find(F==max(F(index_wo_best_model)));

if F(index_best_model) - F(index_second_model) > 3
    winning = 1;
    disp(['Model ', all_model_architecture{index_best_model}, ' is the winning model.'])
else
    winning = 0;
    disp(['No winning model :('])
end

end