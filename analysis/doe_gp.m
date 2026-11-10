% Gaussian-process surrogate of the patterned-wafer DOE (MATLAB R2025b, Statistics and Machine Learning Toolbox).
% Fits log10(t90_worst), edge_lag_pct and utilization on the coded factors, leave-one-out R^2, then finds the
% operating point that minimises the saturating pulse subject to edge lag <= 5 %, and reports it with the GP's
% uncertainty. Also Monte Carlo over the model constants s and N0 (log-normal, factor 1.5 and 1.25) through the
% supply-limit rule to put bands on the baseline and optimum.
root = '/Users/thegreat/Documents/MATLAB/ald-conformality/';
T = readtable([root 'figures/doe_results.csv']);
X = [T.c_x T.c_Q T.c_p T.c_g];
resp = struct('name', {'log10_t90_worst', 'edge_lag_pct', 'utilization'}, 'y', {log10(T.t90_worst), T.edge_lag_pct, T.utilization});
models = cell(1, 3); fid = fopen([root 'figures/doe_gp_report.txt'], 'w');
fprintf(fid, 'MATLAB %s GP surrogates on %d DOE runs\n', version, height(T));
for k = 1:3
    y = resp(k).y; ok = ~isnan(y);
    gp = fitrgp(X(ok, :), y(ok), 'KernelFunction', 'ardsquaredexponential', 'BasisFunction', 'linear', 'Standardize', true);
    % leave-one-out
    n = sum(ok); yhat = zeros(n, 1); Xo = X(ok, :); yo = y(ok);
    for i = 1:n
        idx = true(n, 1); idx(i) = false;
        g = fitrgp(Xo(idx, :), yo(idx), 'KernelFunction', 'ardsquaredexponential', 'BasisFunction', 'linear', 'Standardize', true);
        yhat(i) = predict(g, Xo(i, :));
    end
    r2 = 1 - sum((yo - yhat).^2) / sum((yo - mean(yo)).^2);
    fprintf(fid, '%s: LOO R^2 = %.3f, LOO RMSE = %.4f\n', resp(k).name, r2, sqrt(mean((yo - yhat).^2)));
    models{k} = gp;
end
% grid search over the coded cube for the fastest pulse with edge lag <= 5 %
g = linspace(-1, 1, 21); [A, B, C, D] = ndgrid(g, g, g, g); G = [A(:) B(:) C(:) D(:)];
[t_pred, t_sd] = predict(models{1}, G); lag_pred = predict(models{2}, G); util_pred = predict(models{3}, G);
feas = lag_pred <= 5;
[~, ib] = min(t_pred + 1e6 * ~feas);
[~, i0] = min(sum(abs(G), 2));                      % centre point = baseline-like
dec = @(c) [0.005 * (0.03/0.005).^((c(1)+1)/2) * 100, 150 * 4.^((c(2)+1)/2), 100 * 10.^((c(3)+1)/2), 10 + 7.5 * (c(4)+1)];
fprintf(fid, '\nCentre point: x=%.2f %%, Q=%.0f sccm, p=%.0f Pa, g=%.1f mm -> t90 %.2f s (sd factor %.2f), edge lag %.1f %%, utilization %.0f %%\n', dec(G(i0,:)), 10^t_pred(i0), 10^t_sd(i0), lag_pred(i0), 100*util_pred(i0));
fprintf(fid, 'Fastest with edge lag <= 5 %%: x=%.2f %%, Q=%.0f sccm, p=%.0f Pa, g=%.1f mm -> t90 %.2f s (sd factor %.2f), edge lag %.1f %%, utilization %.0f %%\n', dec(G(ib,:)), 10^t_pred(ib), 10^t_sd(ib), lag_pred(ib), 100*util_pred(ib));
% same precursor supply as the centre (x*Q fixed): best uniformity and utilization
same_supply = abs((G(:,1) + G(:,2)) - (G(i0,1) + G(i0,2))) < 1e-9;   % log-coded x and Q both span factor 6 and 4: not equal; use product in real units
xq = (0.005 * (0.03/0.005).^((G(:,1)+1)/2)) .* (150 * 4.^((G(:,2)+1)/2)); xq0 = xq(i0);
same_supply = abs(log(xq / xq0)) < 0.05;
cand = find(same_supply & feas); [~, j] = min(t_pred(cand)); ib2 = cand(j);
fprintf(fid, 'Same precursor supply as centre, edge lag <= 5 %%, fastest: x=%.2f %%, Q=%.0f sccm, p=%.0f Pa, g=%.1f mm -> t90 %.2f s, edge lag %.1f %%, utilization %.0f %%\n', dec(G(ib2,:)), 10^t_pred(ib2), lag_pred(ib2), 100*util_pred(ib2));
% Monte Carlo on the rule constants for the centre point
rng(1); N = 10000; s = 5e-3 * 1.5.^randn(N,1); N0 = 4e18 * 1.25.^randn(N,1);
tau = 1.2; wall = 1.78e13 * pi * 100e-9 * 100 * 100e-9; Uc0 = 0.0501 * (0.0122 * 300 / (1.380649e-23 * 573.15));
t_rule = tau * (1 + wall) .* N0 / Uc0;
fprintf(fid, '\nMonte Carlo (N0 lognormal x/÷1.25, 10000 draws) on the rule at the centre point: t_pulse = %.1f s, 90 %% band %.1f to %.1f s\n', median(t_rule), prctile(t_rule,5), prctile(t_rule,95));
fclose(fid); type([root 'figures/doe_gp_report.txt']);
