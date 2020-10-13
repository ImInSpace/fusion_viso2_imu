function print_vec(J,name,X,U,do_simplify)
    if nargin<5
        do_simplify=0;
    end

    
    fprintf("VectorXd %s = VectorXd::Zero(%i);\n",name,length(J))
    if do_simplify
        J=simplify(J,100);
        
        if nargin>2 
            syms x(t) u(t)
            J=subs(J,X,sym('x',size(X)));
            J=subs(J,U,sym('u',size(U)));
        end
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
            tline=regexprep(tline,'A0\[(\d*)\]\[0\]',[name '($1)']);
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
        for i=1:length(J)
                if J(i)~=0
                    Ji=sprintf('%s(%i)',name,i-1);
                    print_var(J(i),Ji,X,U);
                end
        end
    end
        fprintf('\n');
end
