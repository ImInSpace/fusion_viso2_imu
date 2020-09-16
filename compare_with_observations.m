addpath(genpath('third-party/yamlmatlab'))
figure('WindowState', 'maximized')
gt_file = 'vehicle_data/GT.csv';
u_file = 'vehicle_data/U.csv';
obs_files = {'vehicle_data/IMU.csv','vehicle_data/VISUAL.csv'};

tc=3;

config=yaml.ReadYaml('config.yaml');

if isfield(config,'GT_from_file')
    config=rmfield(config,'GT_from_file');
end
if isfield(config,'U_from_file')
    config=rmfield(config,'U_from_file');
end
config.GT_to_file=gt_file;
config.U_to_file=u_file;
for i=1:2
    if isfield(config.testCase{tc}.obs{i},'from_file')
        config.testCase{tc}.obs{i}=rmfield(config.testCase{tc}.obs{i},'from_file');
    end
    config.testCase{tc}.obs{i}.to_file=obs_files{i};
end

configObs=config.testCase{tc}.obs;
config.testCase{tc}.obs(:)=[];


yaml.WriteYaml('config2.yaml',config);
subplot(2,2,1);
exec_and_plot();


config=rmfield(config,'GT_to_file');
config=rmfield(config,'U_to_file');
config.GT_from_file=gt_file;
config.U_from_file=u_file;


config.testCase{tc}.obs(1) = configObs(1);
yaml.WriteYaml('config2.yaml',config);
subplot(2,2,3);
exec_and_plot()


config.testCase{tc}.obs(1) = configObs(2);

subplot(2,2,2);

yaml.WriteYaml('config2.yaml',config);
exec_and_plot([2 1 3 4])

config.testCase{tc}.obs(1:2) = configObs(1:2);

for i=1:2
    config.testCase{tc}.obs{i}=rmfield(config.testCase{tc}.obs{i},'to_file');
    config.testCase{tc}.obs{i}.from_file=obs_files{i};
end

yaml.WriteYaml('config2.yaml',config);
subplot(2,2,4);
exec_and_plot()


function exec_and_plot(obs_idx)
    if nargin==0
        obs_idx=1:4;
    end
    system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config2.yaml');
    X=csvread("data.csv");
    tx=X(:,1); ty=X(:,2);
    nObs=X(1,3);
    ox=X(:,4:2:(2+2*nObs)); oy=X(:,5:2:(3+2*nObs));
    kx=X(:,4+2*nObs); ky=X(:,5+2*nObs);
    Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);

    hold on
    scatter(tx(1),ty(1),30,'r');
    legend_text={'start'};
    plot(tx,ty,'r');
    legend_text{end+1}='ground\_truth';
    obs_colors='gcmy';
    for i=1:nObs
        idx=obs_idx(i);
        plot(ox(:,i),oy(:,i),obs_colors(idx))
        legend_text{end+1}=['observation ' num2str(idx)];
    end
    plot(kx,ky,'b')
    legend_text{end+1}='kalman';
    axis equal
    if nObs==1
        obs_legend={'observation'};
    end
    lgd = legend(legend_text);
    lgd.Location='Best';
    
    dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
    dist_error=vecnorm([kx-tx ky-ty],2,2);
    dlm = fitlm(dist_traveled,dist_error,'Intercept',false);
    title({sprintf('Kalman Drift: %.1f%%',dlm.Coefficients.Estimate*100),...
           sprintf('Final kalman uncertainty: %.5f',nthroot(abs(det(Pk(:,:,end))),4))})
    drawnow()
end