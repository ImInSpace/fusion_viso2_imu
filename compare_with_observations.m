addpath(genpath('third-party/yamlmatlab'))
addpath(genpath('matlab_functions'))
figure('WindowState', 'maximized')
gt_file = 'vehicle_data/GT.csv';
u_file = 'vehicle_data/U.csv';
obs_files = {'vehicle_data/IMU.csv','vehicle_data/VISUAL.csv'};

tc=3;
axes=[];

config=yaml.ReadYaml('config.yaml');
configCase=config.testCase{tc};


if isfield(configCase,'GT_from_file')
    configCase=rmfield(configCase,'GT_from_file');
end
if isfield(configCase,'U_from_file')
    configCase=rmfield(configCase,'U_from_file');
end
configCase.GT_to_file=gt_file;
configCase.U_to_file=u_file;
for i=1:2
    if isfield(configCase.obs{i},'from_file')
        configCase.obs{i}=rmfield(configCase.obs{i},'from_file');
    end
    configCase.obs{i}.to_file=obs_files{i};
end

configObs=configCase.obs;
configCase.obs(:)=[];


config.testCase{tc}=configCase;
yaml.WriteYaml('config2.yaml',config);
axes(end+1)=subplot(2,2,1);
exec_and_plot();


configCase=rmfield(configCase,'GT_to_file');
configCase=rmfield(configCase,'U_to_file');
configCase.GT_from_file=gt_file;
configCase.U_from_file=u_file;


configCase.obs(1) = configObs(1);
config.testCase{tc}=configCase;
yaml.WriteYaml('config2.yaml',config);
axes(end+1)=subplot(2,2,3);
exec_and_plot()


configCase.obs(1) = configObs(2);

axes(end+1)=subplot(2,2,2);

config.testCase{tc}=configCase;
yaml.WriteYaml('config2.yaml',config);
exec_and_plot([2 1 3 4])

configCase.obs(1:2) = configObs(1:2);

for i=1:2
    configCase.obs{i}=rmfield(configCase.obs{i},'to_file');
    configCase.obs{i}.from_file=obs_files{i};
end

config.testCase{tc}=configCase;
yaml.WriteYaml('config2.yaml',config);
axes(end+1)=subplot(2,2,4);
exec_and_plot()

linkaxes(axes,'xy')

function exec_and_plot(obs_idx)
    if nargin==0
        obs_idx=1:4;
    end
    exec_and_plot_dets(obs_idx)
end


function exec_and_plot_dets(obs_idx)
    system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config2.yaml');
    X=csvread("data.csv");
    Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);
    dets=arrayfun(@(i)det(Pk(:,:,i)),1:size(Pk,3));
    plot(nthroot(abs(dets),4))
end

function exec_and_plot_locs(obs_idx)
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