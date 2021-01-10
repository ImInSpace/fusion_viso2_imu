vd=[];
id=[];
fd=[];
%close all

dims=4;

%fig=figure('WindowState','maximized');
%tiles=tiledlayout(3,5, 'Padding', 'none', 'TileSpacing', 'compact'); 
is=1;
for i = is
%ax=nexttile;
if dims<6
    ax=figure;
    hold on
end
for old = 0
    if old
base = "\\wsl$\docker-desktop-data\version-pack-data\community\docker\volumes\exoter_logs\_data\processed_old\log-";
    else
base = "\\wsl$\docker-desktop-data\version-pack-data\community\docker\volumes\exoter_logs\_data\processed\log-";
    end
logdir=strcat(base,num2str(i));
if logdir(end)~='\' || logdir(end)~='\'
    logdir=strcat(logdir,'\');
end
addpath(genpath('third-party/yamlmatlab'))
addpath(genpath('matlab_functions'))

matpath=strcat(logdir,'data.mat');
if exist(matpath)
    load(matpath)
else
    TT_VO=read_log(strcat(logdir,'eval\visual_in_wp.txt'));
    TT_IO=read_log(strcat(logdir,'eval\inertial_in_wp.txt'));
    TT_KF=read_log(strcat(logdir,'eval\fusion_in_wp.txt'));
    TT_GT=read_log(strcat(logdir,'gt.txt'));
    TT_VO2=read_log(strcat(logdir,'vo.txt'));
    TT_IO2=read_log(strcat(logdir,'io.txt'));
    TT_KF2=read_log(strcat(logdir,'kf.txt'));
    F_PE=read_log(strcat(logdir,'eval\fusion_pos_err.txt'));
    F_DT=read_log(strcat(logdir,'eval\fusion_dist.txt'));
    V_PE=read_log(strcat(logdir,'eval\visual_pos_err.txt'));
    V_DT=read_log(strcat(logdir,'eval\visual_dist.txt'));
    I_PE=read_log(strcat(logdir,'eval\inertial_pos_err.txt'));
    I_DT=read_log(strcat(logdir,'eval\inertial_dist.txt'));
    V_YE=read_log(strcat(logdir,'eval\visual_yaw_err.txt'));
    I_YE=read_log(strcat(logdir,'eval\inertial_yaw_err.txt'));
    F_YE=read_log(strcat(logdir,'eval\fusion_yaw_err.txt'));
    save(strcat(logdir,'data.mat'),'TT_VO','TT_IO','TT_KF','TT_GT','TT_VO2','TT_IO2','TT_KF2',...
        'F_PE','F_DT','V_PE','V_DT','I_PE','I_DT','I_YE','F_YE','V_YE')
end



TT_GT(abs(TT_GT.euler_angles(:,2))>15*pi/180 | abs(TT_GT.euler_angles(:,3))>15*pi/180,:)=[];
    lims=[0 9];
if dims==7
    fd(i)=mean(seconds(diff(TT_KF2.Time)));
    vd(i)=mean(seconds(diff(TT_VO2.Time)));
    id(i)=mean(seconds(diff(TT_IO2.Time)));
elseif dims==6
    fd(i)=mean(seconds(diff(TT_KF2.Time)));
    vd(i)=mean(seconds(diff(TT_VO2.Time)));
    id(i)=mean(seconds(diff(TT_IO2.Time)));
elseif dims==3
axis equal
    scatter3(TT_GT.position(:,1),TT_GT.position(:,2),TT_GT.position(:,3),1,'filled','r')
    scatter3(TT_IO.position(:,1),TT_IO.position(:,2),TT_IO.position(:,3),1,'filled','g')
    scatter3(TT_VO.position(:,1),TT_VO.position(:,2),TT_VO.position(:,3),1,'filled','c')
    %scatter3(TT_KF.position(:,1),TT_KF.position(:,2),TT_KF.position(:,3),1,'filled','b')
    %xlim(lims)
    %ylim(lims)
    legend gt io vo %kf
elseif dims==2
axis equal
    vars={};
    %varargs={'Marker','.','MarkerSize',5};
    %if old
    %    vars=[vars,'LineStyle','--'];
    %else
    plot(TT_GT.position(2:end,1),TT_GT.position(2:end,2),'r',vars{:})
    plot(TT_IO.position(:,1),TT_IO.position(:,2),'g',vars{:})
    plot(TT_VO.position(:,1),TT_VO.position(:,2),'c',vars{:})
    %end
    plot(TT_KF.position(:,1),TT_KF.position(:,2),'b',vars{:})
    drawnow();
    xlim(lims)
    ylim(lims)
    xlabel x[m]
    ylabel y[m]
    grid on
    lgd=legend({'truth','inertial','visual','kalman'});
    lgd.Location='northwest';
else
    if dims==5
        F_PE=F_YE(1:height(F_DT),:);
        I_PE=I_YE(1:height(I_DT),:);
        V_PE=V_YE(1:height(V_DT),:);
    else
        F_PE=F_PE(1:height(F_DT),:);
        I_PE=I_PE(1:height(I_DT),:);
        V_PE=V_PE(1:height(V_DT),:);
    end
    vars={};
    %varargs={'Marker','.','MarkerSize',5};
    %if old
    %    vars=[vars,'LineStyle','--'];
    %else
    plot(I_DT.Variables,I_PE.Variables,'g',vars{:})
    plot(V_DT.Variables,V_PE.Variables,'c',vars{:})
    %end
    plot(F_DT.Variables,F_PE.Variables,'b',vars{:})
    plot(F_DT.Variables,0.02*F_DT.Variables,'black','LineStyle','--',vars{:})
    
    xlim([0 max([max(F_DT.Variables),max(V_DT.Variables),max(I_DT.Variables)])])
    
    dlm = fitlm(I_DT.Variables,I_PE.Variables,'Intercept',false);
    id(i)=dlm.Coefficients.Estimate*100;
    fprintf('Inertial: %.2f%%\n',dlm.Coefficients.Estimate*100)
    
    dlm = fitlm(V_DT.Variables,V_PE.Variables,'Intercept',false);
    vd(i)=dlm.Coefficients.Estimate*100;
    fprintf('Visual:   %.2f%%\n',dlm.Coefficients.Estimate*100)
    
    
    xlabel("Distance travelled [m]")
    if dims==5
    ylabel("Yaw error [rad]")
    lgd=legend({'inertial','visual','kalman'});
    else
    ylabel("Position error [m]")
    dlm = fitlm(F_DT.Variables,F_PE.Variables,'Intercept',false);
    fd(i)=dlm.Coefficients.Estimate*100;
    fprintf('Fusion:   %.2f%%\n',dlm.Coefficients.Estimate*100)
    lgd=legend({'inertial','visual','fusion','2% drift'});
    
    end
    lgd.Location='northwest';
end
end
if dims<6
    grid on
    %title(strcat('log-',num2str(i)))
    drawnow
    file=strcat('log_',num2str(i),'.pdf');
    %exportgraphics(ax,file)
end
end
if dims==6
id=1/mean(id);vd=1/mean(vd);fd=1/mean(fd);
end
js=1:length(id);
f=js~=7 & js~=10 & id~=0;

if dims==4 || dims==6
   for j=1:length(id)
       ds=[id(j) vd(j) fd(j)];
       if dims==4
            mind=find(min(ds)==ds);
       else
            mind=find(max(ds)==ds);
       end
       hs=["" "" ""];
       hs(mind)="\col";
       fprintf('log-%i & %s\\SI{%.2f}{\\%%} & %s\\SI{%.2f}{\\%%} & %s\\SI{%.2f}{\\%%} \\\\',...
                     j, hs(1), id(j),    hs(2), vd(j),    hs(3), fd(j))
       if j==7 || j==10
           disp("[-1.5ex]")
           disp("\hline\noalign{\vspace{\dimexpr 1.5ex-\doublerulesep}}")
       else
           fprintf("\n")
       end
   end
end
figure;
if length(id)>1
    boxplot([id(f).' vd(f).' fd(f).'],'Labels',{'Inertial','Visual','Filter'})
    ylabel('Drift [%]')
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    [pr,hr]=signrank(vd(f),fd(f),'tail','right');
    disp(['Probablity that fusion is better: ' num2str((1-pr)*100)])
end