syms mu
X=sym('X',[4,1]);
r1=((X(1)-mu)^2+X(2)^2)^(1/2);
r2=((X(1)-mu+1)^2+X(2)^2)^(1/2);
S=.5*(X(1)^2+X(2)^2)+(1-mu)/r1+mu/r2+.5*mu*(1-mu);
Sx=diff(S,X(1));
Sy=diff(S,X(2));
F=[X(3);
    X(4);
    Sx+2*X(4);
    Sy-2*X(3)];
J=jacobian(F,X);
r1_s=sym('r1');
r2_s=sym('r2');

F=subs(F,[r1 r2],[r1_s r2_s])
J=subs(J,[r1 r2],[r1_s r2_s])