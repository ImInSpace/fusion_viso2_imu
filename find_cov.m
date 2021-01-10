eio=[];evo=[];egt=[];
for i = 1:15
    load(strcat('exoter_logs\log-',num2str(i),'\filt.mat'))
    dio=[TT_DIO.velocity TT_DIO.angular_velocity];
    gtio=GT_IO(:,7:12);
    dio=dio(1:size(gtio,1),:);
    neio=dio-gtio;
    eio=[eio;neio(10:end-10,:)];
    
    egt=[egt;diff(GT_IO(:,7:12))./repmat(seconds(diff(TT_DIO.Time(1:size(gtio,1)))),1,6)];
    
    gtvo=GT_VO(:,1:6);
    dgtvo=zeros(size(gtvo,1)-1,6);
    for j=2:size(gtvo,1)
        pos_s=gtvo(j-1,1:3).';
        pos=gtvo(j,1:3).';
        eul_s=gtvo(j-1,4:6).';
        eul=gtvo(j,4:6).';
        R=my_eul2rotm(eul_s);
        euld2lw=euld2localw(eul_s);
        dgtvo(j-1,1:3)=(R.'*(pos-pos_s)).';
        dgtvo(j-1,4:6)=(euld2lw*angdiff(eul_s,eul)).';
    end
    dvo=[TT_DVO.position flip(TT_DVO.euler_angles,2)];
    dvo=dvo(2:end,:);
    dvo=dvo(1:size(dgtvo,1),:);
    nevo=dvo-dgtvo;
    evo=[evo;nevo(4:end-4,:)];
end
evo = rmoutliers(evo);
eio = rmoutliers(eio);
egt = rmoutliers(egt);
Cio=cov(eio);
stdio=std(eio);
Cvo=cov(evo);
stdvo=std(evo);
stdgt=std(egt);