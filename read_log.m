function read_log(logfile)
    opts=detectImportOptions(logfile); %detect automatic import options
    %If you look at the variablenames here you see that they are shifted
    %from their corresponding data This is because sourceFrame is a string
    %with spaces and it get detected as more columns
    
    %detect where sourceframe is
    sourceframe_col=find(contains(opts.VariableNames,'sourceFrame'));
    %detect how many extra columns there are
    num_extra_cols=sum(startsWith(opts.VariableNames,'Var'));
    
    % shift everything to their correct spot
    opts.VariableNames(sourceframe_col+num_extra_cols:end)=opts.VariableNames(sourceframe_col:end-num_extra_cols);
    % put new names to the extra columns of sourceframe
    opts.VariableNames(sourceframe_col:sourceframe_col+num_extra_cols-1)=arrayfun(@(i)[opts.VariableNames{sourceframe_col} num2str(i)],1:num_extra_cols,'UniformOutput',false);
    
    %remove x_vicon_pose_samples_ from the start of all variable names
    opts.VariableNames=cellfun(@(name)name(22:end),opts.VariableNames,'UniformOutput',false);
    %remove underscore from ending of some variables
    end_with_subs=endsWith(opts.VariableNames,'_');
    opts.VariableNames(end_with_subs)=cellfun(@(name)name(1:end-1),opts.VariableNames(end_with_subs),'UniformOutput',false);
    %remove __1 from ending of some variables
    end_with_1=endsWith(opts.VariableNames,'__1');
    opts.VariableNames(end_with_1)=cellfun(@(name)name(1:end-3),opts.VariableNames(end_with_1),'UniformOutput',false);
    
    %read the file
    T=readtable(logfile,opts);
    %remove nan rows
    T(isnan(T.cov_position_data_1),:)=[]; 
    
    %convert to timetable
    TT=timetable(datetime(T.time_microseconds*10^-6,'ConvertFrom','posixtime'));
    
    %add time elapsed
    TT=addvars(TT,seconds(TT.Time-TT.Time(1)),'NewVariableNames','elapsed_seconds');
    
    %add a variable with quaternion
    TT=addvars(TT,quaternion(T.orientation_re,T.orientation_im_0,T.orientation_im_1,T.orientation_im_2),'NewVariableNames','quaternion');
    %add a variable with euler angles
    TT=addvars(TT,quat2eul(TT.quaternion),'NewVariableNames','euler_angles');
    %format position
    TT=addvars(TT,[T.position_data_0 T.position_data_1 T.position_data_2],'NewVariableNames','position');
    
    %TT(abs(TT.euler_angles(:,2))>15*pi/180 | abs(TT.euler_angles(:,3))>15*pi/180,:)=[];
  
    
    assignin('base','TT',TT)
    
    figure()
    hold on 
    
    axis equal
    
    p=scatter3(TT.position(:,1),TT.position(:,2),TT.position(:,3),2,'filled');
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('psi',TT.euler_angles(:,1));
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('theta',TT.euler_angles(:,2));
    p.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('phi',TT.euler_angles(:,3));
    
    %plot3(x,y,z);
    
    car1=stlread('car_1.stl');
    V=car1.Points/3;
    F=car1.ConnectivityList;
    V(:,3)=-V(:,3);
    V(:,1)=-V(:,1);
    pat=patch('Vertices',V,'Faces',F,'FaceColor','r');
    tic
    for i=ceil(linspace(1,size(TT,1),200))
        while toc<TT.elapsed_seconds(i)/100
        end
        R=eul2rotm(TT.euler_angles(i,:));
        %R=rotmat(T.quaternion(i),'point');
        P=TT.position(i,:);
        pat.set('Vertices',P+V*R.')
        drawnow
    end
    
    
    csvwrite(regexprep(logfile,'\.txt','.csv'),[TT.position TT.euler_angles]);
    csvwrite(regexprep(logfile,'\.txt','_time.csv'),TT.elapsed_seconds);
end