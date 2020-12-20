fid = fopen('sa.txt');
tline = fgetl(fid);
vn=split(tline,' ');
len_prefix=length(vn{1})-16;
vn=cellfun(@(name)name(len_prefix:end),vn,'UniformOutput',false);
    
%disp(length(a))
while ischar(tline)
    data=split(tline,' ');
    data=cellfun(@(d)help(str2double(d)),data,'UniformOutput',false);
    ndata={};
    for i=1:length(data)
        if ischar(data{i})
            if ischar(ndata{end})
                ndata{end}=[ndata{end} data{i}];
            else
                ndata{end+1}=data{i};
            end
        else
            ndata{end+1}=data{i};
        end
    end
    
    %disp(length(a))
    tline = fgetl(fid);
end

function res=help(a)
    if a>=56 && a<=172
        res=char(a);
    else
        res=a;
    end
end