function print_vec(F,name,X,U)
    syms x(t) u(t)
    F=subs(F,X,x(0:length(X)-1).');
    F=subs(F,U,u(0:length(U)-1).');
    fprintf("VectorXd %s(%i);\n",name,length(F))
    fprintf("%s << ",name);
    for i=1:length(F)
        Fi=ccode(F(i));
        Fi=regexprep(Fi,'  t0 = ','');
        Fi=regexprep(Fi,';','');
        Fi=regexprep(Fi,'.0)',')');
        if(i==length(F))
            fprintf('%s;\n\n',Fi)
        else
            fprintf('%s,\n',Fi)
        end
    end
end