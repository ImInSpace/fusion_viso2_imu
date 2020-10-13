function read_log(logfile)
    opts=detectImportOptions(logfile);
    sourceframe_col=find(contains(opts.VariableNames,'sourceFrame'));
    num_extra_cols=sum(startsWith(opts.VariableNames,'Var'));
    
    opts.VariableNames(sourceframe_col+num_extra_cols:end)=opts.VariableNames(sourceframe_col:end-num_extra_cols);
    opts.VariableNames(sourceframe_col:sourceframe_col+num_extra_cols-1)=arrayfun(@(i)[opts.VariableNames{sourceframe_col} num2str(i)],1:num_extra_cols,'UniformOutput',false);
    opts.VariableNames=cellfun(@(name)name(22:end),opts.VariableNames,'UniformOutput',false);
    end_with_subs=endsWith(opts.VariableNames,'_');
    opts.VariableNames(end_with_subs)=cellfun(@(name)name(1:end-1),opts.VariableNames(end_with_subs),'UniformOutput',false);
    end_with_1=endsWith(opts.VariableNames,'__1');
    opts.VariableNames(end_with_1)=cellfun(@(name)name(1:end-3),opts.VariableNames(end_with_1),'UniformOutput',false);
    
    T=readtable(logfile,opts);
    T(isnan(T.cov_position_data_1),:)=[];
    T=addvars(T,quaternion(T.orientation_re,T.orientation_im_0,T.orientation_im_1,T.orientation_im_2),'NewVariableNames','quaternion');
    
    figure()
    hold on 
    
    axis equal
    
    p=scatter3(T.position_data_0,T.position_data_1,T.position_data_2,2,'filled');
    pl=plot3(T.position_data_0,T.position_data_1,T.position_data_2);
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cov_x',T.cov_position_data_0);
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cov_y',T.cov_position_data_4);
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('cov_z',T.cov_position_data_8);
    
    x=filter(ones(1,5)/5,1,T.position_data_0);
    y=filter(ones(1,5)/5,1,T.position_data_1);
    z=filter(ones(1,5)/5,1,T.position_data_2);
    %plot3(x,y,z);
    
    car1=stlread('car_1.stl');
    V=car1.Points/3;
    F=car1.ConnectivityList;
    V(:,3)=-V(:,3);
    V(:,1)=-V(:,1);
    pat=patch('Vertices',V,'Faces',F,'FaceColor','r');
    t0=T.time_microseconds(1);
    tic
    for i=ceil(linspace(1,size(T,1),2000))
        while toc<(T.time_microseconds(i)-t0)*10^-6/10
        end
        eul=quat2eul(T.quaternion(i));
        R=eul2rotm(eul);
        %R=rotmat(T.quaternion(i),'point');
        P=[T.position_data_0(i) T.position_data_1(i) T.position_data_2(i)];
        pat.set('Vertices',P+V*R.')
        drawnow
    end
    
    TT=timetable(datetime(T.time_microseconds*10^-6,'ConvertFrom','posixtime'));
    assignin('base','TT',TT)
    assignin('base','T',T)
    assignin('base','opts',opts)
    assignin('base','logfile',logfile)
end