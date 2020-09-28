
addpath(genpath('matlab_functions'))
figure()
system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config.yaml');
plot_data
lgd_mix=cellfun(@(l)['vel\_' l], legend_text, 'UniformOutput', false);

dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
dist_error=vecnorm([kx-tx ky-ty],2,2);
dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
fprintf('Vel Kalman Drift: %.1f\n',dist_traveled(end)*dlm.Coefficients.Estimate)
fprintf('Vel Final kalman uncertainty: %.5f\n',nthroot(abs(det(Pk(:,:,end))),4))


system('cmake-build-debug\test_fusion_viso2_imu.exe 4 0 config.yaml');
plot_data
lgd_mix=[lgd_mix cellfun(@(l)['clo\_' l], legend_text, 'UniformOutput', false)];
legend(lgd_mix)

hLines = findobj(gca,'Type','line');
for i=1:length(hLines)
    line_i=hLines(i);
    if startsWith(line_i.DisplayName,'clo')
        line_i.LineStyle='--';
        if ~endsWith(line_i.DisplayName,'kalman')
            delete(line_i)
        end
    end
end


hScatters = findobj(gca,'Type','Scatter');
for i=1:length(hScatters)
    scatter_i=hScatters(i);
    if strcmp(scatter_i.DisplayName,'clo\_start')
        delete(scatter_i)
    end
end

dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
dist_error=vecnorm([kx-tx ky-ty],2,2);
dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
fprintf('Clone Kalman Drift: %.1f\n',dist_traveled(end)*dlm.Coefficients.Estimate)
fprintf('Clone Final kalman uncertainty: %.5f\n',nthroot(abs(det(Pk(:,:,end))),4))