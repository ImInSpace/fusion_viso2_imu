function read_log(logfile)
    opts=detectImportOptions(logfile); %detect automatic import options
    %If you look at the variablenames here you see that they are shifted
    %from their corresponding data This is because sourceFrame is a string
    %with spaces and it get detected as more columns
    
    %detect where sourceframe is
    sourceframe_col=find(contains(opts.VariableNames,'sourceFrame'));
    %detect how many extra columns there are
    num_extra_cols=sum(startsWith(opts.VariableNames,'Var'));
    
    if (num_extra_cols>0)
        % shift everything to their correct spot
        opts.VariableNames(sourceframe_col+num_extra_cols:end)=opts.VariableNames(sourceframe_col:end-num_extra_cols);
        % put new names to the extra columns of sourceframe
        opts.VariableNames(sourceframe_col:sourceframe_col+num_extra_cols-1)=arrayfun(@(i)[opts.VariableNames{sourceframe_col} num2str(i)],1:num_extra_cols,'UniformOutput',false);
    end
    Tprev=preview(logfile,opts);
    if iscell(Tprev(1,end-1).Variables)
        opts.VariableNames(sourceframe_col:sourceframe_col+1)=[];
    end
    
    %remove x_vicon_pose_samples_ from the start of all variable names
    len_prefix=length(opts.VariableNames{1})-16;
    opts.VariableNames=cellfun(@(name)name(len_prefix:end),opts.VariableNames,'UniformOutput',false);
    %remove underscore from ending of some variables
    end_with_subs=endsWith(opts.VariableNames,'_');
    opts.VariableNames(end_with_subs)=cellfun(@(name)name(1:end-1),opts.VariableNames(end_with_subs),'UniformOutput',false);
    %remove __1 from ending of some variables
    end_with_1=endsWith(opts.VariableNames,'__1');
    opts.VariableNames(end_with_1)=cellfun(@(name)name(1:end-3),opts.VariableNames(end_with_1),'UniformOutput',false);
    
    %read the file
    T=readtable(logfile,opts);
    %remove nan rows
    T(isnan(T.position_data_1),:)=[]; 
    
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
    
    if any(~isnan(T.velocity_data_0))
        TT=addvars(TT,[T.velocity_data_0 T.velocity_data_1 T.velocity_data_2],'NewVariableNames','velocity');
    end
    
    if any(~isnan(T.angular_velocity_data_0))
        TT=addvars(TT,[T.angular_velocity_data_0 T.angular_velocity_data_1 T.angular_velocity_data_2],'NewVariableNames','angular_velocity');
    end
    
    %remove poses where pitch or roll is greater than 20 degrees (vicon
    %errors)
    TT(abs(TT.euler_angles(:,2))>20*pi/180 | abs(TT.euler_angles(:,3))>20*pi/180,:)=[];
    
    
    assignin('base','TT',TT)
    assignin('base','T',T)
    
    return
    figure()
    view(3)
    hold on 
    
    axis equal
    
    pl=scatter3(TT.position(:,1),TT.position(:,2),TT.position(:,3),2,'filled');
    pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('psi',TT.euler_angles(:,1)*180/pi);
    pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('theta',TT.euler_angles(:,2)*180/pi);
    pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('phi',TT.euler_angles(:,3)*180/pi);
    
    %plot3(x,y,z);
    
    [car,V]=plot_car();
    
    tic
    for i=ceil(linspace(1,size(TT,1),200))
        while toc<TT.elapsed_seconds(i)/100
        end
        R=my_eul2rotm(TT.euler_angles(i,:));
        %R=rotmat(T.quaternion(i),'point');
        p=TT.position(i,:).';
        car.set('Vertices',(p+R*V')')
        drawnow
    end
    
    
    csvwrite(regexprep(logfile,'\.txt','.csv'),[TT.position TT.euler_angles]);
    csvwrite(regexprep(logfile,'\.txt','_time.csv'),TT.elapsed_seconds);
    
    
    assignin('base','pl',pl)
end