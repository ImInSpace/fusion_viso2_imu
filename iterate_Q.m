close all
N=10;
Qms=0.05*10.^linspace(-1,1.5,N);
kerror=zeros(1,N);

figure()
hold on
for i=1:N
    disp(i);
    addpath(genpath('third-party/yamlmatlab'))
    yamls=yaml.ReadYaml('config.yaml');
    yamls{3}.Qm=Qms(i);
    yaml.WriteYaml('config2.yaml',yamls);
    system('cmake-build-debug\test_fusion_viso2_imu.exe 3 0 config2.yaml');
    X=csvread("data.csv");
    tx=X(:,1); ty=X(:,2);
    nObs=X(1,3);
    ox=X(:,4:2:(2+2*nObs)); oy=X(:,5:2:(3+2*nObs));
    kx=X(:,4+2*nObs); ky=X(:,5+2*nObs);
    Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);

    hold on
    scatter(tx(1),ty(1),30,'r');
    plot(tx,ty,'r');
    obs_colors='gcmy';
    for iObs=1:nObs
        plot(ox(:,iObs),oy(:,iObs),obs_colors(iObs))
    end
    plot(kx,ky,'b')
    drawnow
    Pk=reshape(X(:,end-3:end).',[2 2 size(X,1)]);

    dist_traveled=[0;cumsum(sqrt(diff(tx).^2+diff(ty).^2))];
    dist_oerror=vecnorm([ox(:,1)-tx oy(:,1)-ty],2,2);
    dist_kerror=vecnorm([kx-tx ky-ty],2,2);
    dlmo = fitlm(dist_traveled,dist_oerror,'Intercept',false);
    dlmk = fitlm(dist_traveled,dist_kerror,'Intercept',false);
    oe=dlmo.Coefficients.Estimate*100;
    ke=dlmk.Coefficients.Estimate*100;
    fprintf('%.1f%%\n',ke/oe*100);
    kerror(i)=ke/oe*100;
end