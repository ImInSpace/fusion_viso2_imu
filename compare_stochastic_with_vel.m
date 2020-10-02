
addpath(genpath('matlab_functions'))
addpath(genpath('third-party/yamlmatlab'))
figure()


config=yaml.ReadYaml('config.yaml');
vel = config.testCase{3};
clo = config.testCase{4};

if ~isfield(vel,'GT_to_file')
    vel.GT_to_file='vehicle_data/GT.csv';
end
if ~isfield(vel,'U_to_file')
    vel.U_to_file='vehicle_data/U.csv';
end
if ~isfield(vel.obs{1},'to_file')
    vel.obs{1}.to_file='vehicle_data/IMU.csv';
end
if ~isfield(vel.obs{2},'to_file')
    vel.obs{2}.to_file='vehicle_data/IMU.csv';
end
if isfield(vel,'GT_from_file')
    vel=rmfield(vel,'GT_from_file');
end
if isfield(vel,'U_from_file')
    vel=rmfield(vel,'U_from_file');
end
if isfield(vel.obs{1},'from_file')
    vel.obs{1}=rmfield(vel.obs{1},'from_file');
end
if isfield(vel.obs{2},'from_file')
    vel.obs{2}=rmfield(vel.obs{2},'from_file');
end

clo.GT_from_file=vel.GT_to_file;
if isfield(clo,'GT_to_file')
    clo=rmfield(clo,'GT_to_file');
end
clo.U_from_file=vel.obs{1}.to_file;
if isfield(clo,'U_to_file')
    clo=rmfield(clo,'U_to_file');
end
clo.obs{1}.from_file=vel.obs{2}.to_file;
if isfield(clo.obs{1},'to_file')
    clo.obs{1}=rmfield(clo.obs{1},'to_file');
end

clo.dt=vel.dt;
clo.T=vel.T;
clo.x0=vel.x0;
clo.Qd=vel.Qd;
clo.Qm=vel.Qm;
clo.Cd=vel.obs{1}.Rd;
clo.Cm=vel.obs{1}.Rm;
clo.obs{1}.every_X=vel.obs{2}.every_X;
clo.obs{1}.Rd=vel.obs{2}.Rd;
clo.obs{1}.Rm=vel.obs{2}.Rm*vel.dt*vel.obs{2}.every_X;


config.testCase{3} = vel;
config.testCase{4} = clo;


yaml.WriteYaml('config2.yaml',config);

system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config2.yaml');
plot_data
lgd_mix=legend_text;
filter=contains(legend_text,'kalman');
lgd_mix(filter)=cellfun(@(l)['vel\_' l], legend_text(filter), 'UniformOutput', false);

dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
dist_error=vecnorm([kx-tx ky-ty],2,2);
dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
fprintf('Vel Kalman Drift: %.1f\n',dist_traveled(end)*dlm.Coefficients.Estimate)
fprintf('Vel Final kalman uncertainty: %.5f\n',nthroot(abs(det(Pk(:,:,end))),4))


system('cmake-build-debug\test_fusion_viso2_imu.exe 4 0 config2.yaml');
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