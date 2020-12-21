for i=1:15
    logdir=strcat('exoter_logs\log-',num2str(i));
    if logdir(end)~='\' || logdir(end)~='\'
        logdir=strcat(logdir,'\');
    end
    addpath(genpath('third-party/yamlmatlab'))
    
    TT_GT=read_log(strcat(logdir,'gt.txt'));
    TT_DVO=read_log(strcat(logdir,'dvo.txt'));
    TT_DIO=read_log(strcat(logdir,'dio.txt'));
    
    
    config=yaml.ReadYaml('config.yaml');
    configCase=config.testCase{5};
    configCase.GT_from_file=strcat(logdir,'gt.csv');
    configCase.obs{1}.from_file=strcat(logdir,'gt.csv');
    configCase.Time_from_file=strcat(logdir,'gt_time.csv');
    
    configCase.OTime_from_file=strcat(logdir,'dvo_time.csv');
    configCase.Kalman_to_file=strcat(logdir,'gt_vo.csv');
    config.testCase{5}=configCase;
    yaml.WriteYaml('config2.yaml',config);
    system('cmake-build-debug\test_fusion_viso2_imu.exe 5 0 config2.yaml');
    GT_VO=csvread(strcat(logdir,'gt_vo.csv'));
    
    configCase.OTime_from_file=strcat(logdir,'dio_time.csv');
    configCase.Kalman_to_file=strcat(logdir,'gt_io.csv');
    config.testCase{5}=configCase;
    yaml.WriteYaml('config2.yaml',config);
    system('cmake-build-debug\test_fusion_viso2_imu.exe 5 0 config2.yaml');
    GT_IO=csvread(strcat(logdir,'gt_io.csv'));
    save(strcat(logdir,'filt.mat'),'GT_IO','GT_VO','TT_DIO','TT_GT')
    delete(strcat(logdir,'*.csv'))
end