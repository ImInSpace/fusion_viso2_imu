function out=read_log(logfile, do_plot)
    if nargin<2
        do_plot=0;
    end
    if ischar(do_plot)
        do_plot=str2double(do_plot);
    end
    %disp(do_plot)
    plot_car=0;
    %Count number of header lines
    fid = fopen(logfile);
    tline = fgetl(fid);
    headerLines=0;
    while ischar(tline) && startsWith(tline,'pocolog')
        headerLines=headerLines+1;
        tline = fgetl(fid);
    end
    fclose(fid);

    opts=detectImportOptions(logfile,'NumHeaderLines',headerLines); %detect automatic import options
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
    if size(Tprev,2)==1
        out=readtable(logfile,opts);
        return
    end
    if iscell(Tprev(1,end-1).Variables)
        opts.VariableNames(sourceframe_col:sourceframe_col+1)=[];
    end
    
    %remove x_vicon_pose_samples_ from the start of all variable names
    %disp(opts.VariableNames{1})
    len_prefix=length(opts.VariableNames{1})-16;
    opts.VariableNames=cellfun(@(name)name(len_prefix:end),opts.VariableNames,'UniformOutput',false);
    %disp(opts.VariableNames{1})
    %remove underscore from ending of some variables
    end_with_subs=endsWith(opts.VariableNames,'_');
    opts.VariableNames(end_with_subs)=cellfun(@(name)name(1:end-1),opts.VariableNames(end_with_subs),'UniformOutput',false);
    %remove __1 from ending of some variables
    end_with_1=endsWith(opts.VariableNames,'__1');
    opts.VariableNames(end_with_1)=cellfun(@(name)name(1:end-3),opts.VariableNames(end_with_1),'UniformOutput',false);
    opts.VariableNames=replace(opts.VariableNames,"__","_");
    %add a that is missing in some cases
    opts.VariableNames=replace(opts.VariableNames,"_dat_","_data_");
    opts.VariableNames=replace(opts.VariableNames,"_velo_","_velocity_");
    end_with_d=endsWith(opts.VariableNames,'_d');
    opts.VariableNames(end_with_d)=replace(opts.VariableNames(end_with_d),"_d","_data");
    end_with_da=endsWith(opts.VariableNames,'_da');
    opts.VariableNames(end_with_da)=replace(opts.VariableNames(end_with_da),"_da","_data");
    end_with_dat=endsWith(opts.VariableNames,'_dat');
    opts.VariableNames(end_with_dat)=replace(opts.VariableNames(end_with_dat),"_dat","_data");
    %add _0 that is missing in some cases
    is_ang_vel_without_0=strcmp(opts.VariableNames,'angular_velocity_data');
    opts.VariableNames(is_ang_vel_without_0)={'angular_velocity_data_0'};
    
    %read the file
    T=readtable(logfile,opts);
    %remove nan rows
    T(isnan(T.position_data_1),:)=[]; 
    
    %convert to timetable
    TT=timetable(datetime(T.time_microseconds*10^-6,'ConvertFrom','posixtime'));
    
    %add time elapsed
    first_t=TT.Time(1);
    first_t=datetime(first_t.Year,first_t.Month,first_t.Day,7,0,0); %save time from 7 am that day
    TT=addvars(TT,seconds(TT.Time-first_t),'NewVariableNames','elapsed_seconds');
    
    %add a variable with quaternion
    TT=addvars(TT,quaternion(T.orientation_re,T.orientation_im_0,T.orientation_im_1,T.orientation_im_2),'NewVariableNames','quaternion');
    %add a variable with euler angles
    TT=addvars(TT,quat2eul(TT.quaternion),'NewVariableNames','euler_angles');
    %format position
    TT=addvars(TT,[T.position_data_0 T.position_data_1 T.position_data_2],'NewVariableNames','position');
    
    velocity_exists = any(~isnan(T.velocity_data_0));
    if velocity_exists
        TT=addvars(TT,[T.velocity_data_0 T.velocity_data_1 T.velocity_data_2],'NewVariableNames','velocity');
    end
    angular_velocity_exists = any(~isnan(T.angular_velocity_data_0));
    if angular_velocity_exists
        TT=addvars(TT,[T.angular_velocity_data_0 T.angular_velocity_data_1 T.angular_velocity_data_2],'NewVariableNames','angular_velocity');
    end
    
    
    %remove poses where pitch or roll is greater than 20 degrees (vicon
    %errors)
    %TT(abs(TT.euler_angles(:,2))>15*pi/180 | abs(TT.euler_angles(:,3))>15*pi/180,:)=[];
    
    %disp([min(TT.elapsed_seconds)-max(TT.elapsed_seconds)])
    if nargout>0
        out=TT;
        return
    end
    assignin('base','TT',TT)
    assignin('base','T',T)
    

    dlmwrite(regexprep(logfile,'\.txt','.csv'),[TT.position TT.euler_angles],'precision',10);
    dlmwrite(regexprep(logfile,'\.txt','_time.csv'),TT.elapsed_seconds,'precision','%.9f');
    if velocity_exists && angular_velocity_exists  
        dlmwrite(regexprep(logfile,'\.txt','_vel.csv'),[TT.velocity TT.angular_velocity],'precision',10);
    end
    if do_plot==1 %plot data directly from position
    
        %figure()
        view(3)
        hold on 

        axis equal
        pl=scatter3(TT.position(:,1),TT.position(:,2),TT.position(:,3),2,'filled');
        pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('psi',TT.euler_angles(:,1)*180/pi);
        pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('theta',TT.euler_angles(:,2)*180/pi);
        pl.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('phi',TT.euler_angles(:,3)*180/pi);

        %plot3(x,y,z);

        if plot_car
        [car,V]=plot_car();
        end

        tic
        for i=ceil(linspace(1,size(TT,1),200))
            if plot_car
                while toc<TT.elapsed_seconds(i)/100
                end
            end
            R=my_eul2rotm(TT.euler_angles(i,:));
            %R=rotmat(T.quaternion(i),'point');
            p=TT.position(i,:).';
            if plot_car
                car.set('Vertices',(p+R*V')')
            end
            drawnow
        end

        assignin('base','pl',pl)
    elseif do_plot==2 %delta pose
        
        figure()
        view(3)
        hold on 

        axis equal
        
        
        [car,V]=plot_car();

        p=[0 0 0].';
        ps=p;
        pl=plot3(ps(1,:),ps(2,:),ps(3,:));
        eul=[0 0 0].';
        euls=eul;
        t=-5;
        for it=1:size(TT,1)
            new_t=TT.elapsed_seconds(it);
            dt=new_t-t;
            t=new_t;
            lw2euld=localw2euld(eul);
            v =(TT.position(it,:).'/dt);
            w = inv(localw2euld([0 0 0]))*TT.euler_angles(it,:).'/dt;
            R=my_eul2rotm(eul);
            eul=eul+(lw2euld*w)*dt;
            p=p+R*v*dt;
            euls(:,end+1)=eul;
            ps(:,end+1)=p;
            car.set('Vertices',(p+R*V.').')
            pl.set('XData',ps(1,:),'YData',ps(2,:),'ZData',ps(3,:));
            
        end
        
        assignin('base','ps',ps)
        assignin('base','euls',euls)
    end
end