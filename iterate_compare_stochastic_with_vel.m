addpath(genpath('third-party/yamlmatlab'))
copyfile('config.yaml','config_orig.yaml')
wb=waitbar(0,'Iterating');
try
    config=yaml.ReadYaml('config.yaml');
    close all
    N=2;
    R1s=10.^linspace(log10(0.002),log10(0.2),N);
    R2s=10.^linspace(log10(0.002),log10(0.2),N);
    [R1s,R2s]=meshgrid(R1s,R2s);
    R1s=R1s(:);
    R2s=R2s(:);
    imu_drifts=zeros(N*N,1);
    vis_drifts=zeros(N*N,1);
    vel_kdrifts=zeros(N*N,1);
    vel_kuncerts=zeros(N*N,1);
    clo_kdrifts=zeros(N*N,1);
    clo_kuncerts=zeros(N*N,1);    
    mkdir iterate_data
    for iter=1:(N*N)
        config.testCase{3}.obs{1}.Rm=R1s(iter);
        config.testCase{3}.obs{2}.Rm=R2s(iter);
        yaml.WriteYaml('config.yaml',config);
        try
        compare_stochastic_with_vel;
        

        imu_drifts(iter)=imu_drift;
        vis_drifts(iter)=vis_drift;
        vel_kdrifts(iter)=vel_kalman_drift;
        vel_kuncerts(iter)=vel_kalman_uncertainty;
        clo_kdrifts(iter)=clo_kalman_drift;
        clo_kuncerts(iter)=clo_kalman_uncertainty;  
        %exportgraphics(fig,['iterate_data/' num2str(i) '.pdf'])
        close(fig);
        waitbar(iter/(N*N),wb);
        catch e
            warning('Error catched:')
            warning(e.message)
            imu_drifts(iter)=nan;
            vis_drifts(iter)=nan;
            vel_kdrifts(iter)=nan;
            vel_kuncerts(iter)=nan;
            clo_kdrifts(iter)=nan;
            clo_kuncerts(iter)=nan; 
        end
    end

catch e
    copyfile('config_orig.yaml','config.yaml')
    close(wb);
    rethrow(e)
end
