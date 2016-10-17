function MD = compute_SMF_distance(smf_real,smf_est, dist)

idx_mode_real = find(smf_real==max(smf_real));
idx_mode_est = find(smf_est==max(smf_est));

d = [];
for i=1:length(idx_mode_real),
    for j=1:length(idx_mode_est),
        d = [d dist(idx_mode_real(i), idx_mode_est(j))];
    end
end

MD = mean(d);
