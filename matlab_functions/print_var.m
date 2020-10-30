function print_var(var,name,X,U)
    if nargin>2
        syms x(t) u(t)
        var=subs(var,X,x(0:length(X)-1).');
        var=subs(var,U,u(0:length(U)-1).');
    end
    
    w=warning('off','all');
    var=ccode(var);
    warning(w);
    var=regexprep(var,'  t0 = ','');
    var=regexprep(var,';','');
    var=regexprep(var,'.0)',')');
    var=regexprep(var,'(\d).(\d)E\+1','$1$2');
    
    fprintf('%s=%s;\n',name,var)
end

