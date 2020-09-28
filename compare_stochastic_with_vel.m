
addpath(genpath('matlab_functions'))
figure()
system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config.yaml');
plot_data
lgd_mix=legend_text;
filter=contains(legend_text,'kalman');
lgd_mix(filter)=cellfun(@(l)['vel\_' l], legend_text(filter), 'UniformOutput', false);

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
    hLine=hLines(i);
    if startsWith(hLine.DisplayName,'clo')
        hLine.LineStyle='--';
        if ~endsWith(hLine.DisplayName,'kalman')
            delete(hLine)
        end
    end
end


hScatters = findobj(gca,'Type','Scatter');
for i=1:length(hScatters)
    hScatter=hScatters(i);
    if strcmp(hScatter.DisplayName,'clo\_start')
        delete(hScatter)
    end
end

hPolys=findobj(gca,'Type','Polygon');

for i=1:length(hPolys)
    hPoly=hPolys(i);
    if startsWith(hPoly.DisplayName,'vel')
        hPoly.LineStyle='-';
    end
end

dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
dist_error=vecnorm([kx-tx ky-ty],2,2);
dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
fprintf('Clone Kalman Drift: %.1f\n',dist_traveled(end)*dlm.Coefficients.Estimate)
fprintf('Clone Final kalman uncertainty: %.5f\n',nthroot(abs(det(Pk(:,:,end))),4))