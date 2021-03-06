function print_mat(J,name,X,U,do_simplify)
    if nargin<5
        do_simplify=0;
    end
    stack=dbstack('-completenames', 1);
    if ~isempty(stack)
        fprintf("  // From %s \n",erase(stack(end).file,pwd))
    end
    fprintf("  MatrixXd %s = MatrixXd::Zero(%i,%i);\n",name,size(J,1),size(J,2))
    if do_simplify
        if nargin>2 
            syms x(t) u(t)
            J=subs(J,X,sym('x',size(X)));
            J=subs(J,U,sym('u',size(U)));
        end
        J=rewrite(simplify(J,100),'cos');
        ccode(J,'File','tempF.cc')
        
        fid = fopen('tempF.cc');
        tline = fgetl(fid);
        while ischar(tline)
            for i=length(X):-1:1
                tline=strrep(tline,['x' num2str(i)],['x(' num2str(i-1) ')']);
            end
            for i=length(U):-1:1
                tline=strrep(tline,['u' num2str(i)],['u(' num2str(i-1) ')']);
            end
            tline=regexprep(tline,'(t\d*) =','double $1 =');
            tline=regexprep(tline,'A0\[(\d*)\]\[(\d*)\]',[name '($1,$2)']);
            disp(tline)
            tline = fgetl(fid);
        end
        fclose(fid);
        delete('tempF.cc')
    else
        if nargin>2 
            syms x(t) u(t)
            J=subs(J,X,x(0:length(X)-1).');
            J=subs(J,U,u(0:length(U)-1).');
        end
        for i=1:size(J,1)
            for j=1:size(J,2)
                if J(i,j)~=0
                    Jij=sprintf('%s(%i,%i)',name,i-1,j-1);
                    print_var(J(i,j),Jij);
                end
            end
        end
    end
    fprintf('  return %s;\n',name)
    fprintf('\n');
end