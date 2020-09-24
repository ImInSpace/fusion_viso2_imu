function print_mat(J,name,X,U)
    syms x(t) u(t)
    J=subs(J,X,x(0:length(X)-1).');
    J=subs(J,U,u(0:length(U)-1).');
    fprintf("MatrixXd %s = MatrixXd::Zero(%i, %i);\n",name,size(J,1),size(J,2))
    for i=1:size(J,1)
        for j=1:size(J,2)
            if J(i,j)~=0
                Jij=ccode(J(i,j));
                Jij=regexprep(Jij,'  t0 = ','');
                Jij=regexprep(Jij,';','');
                Jij=regexprep(Jij,'.0)',')');
                fprintf('%s(%i,%i)=%s;\n',name,i-1,j-1,Jij)
            end
        end
    end
    fprintf('\n');
end