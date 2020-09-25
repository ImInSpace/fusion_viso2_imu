
N=length(X);

J=jacobian(F,X);
H=jacobian(h,X);

%% Discretize
syms dt
B=jacobian(F,U);
%AdBd=expm([J B;zeros(size(B)).' zeros(size(B,2))]*dt);
%A=AdBd(1:size(J,1),1:size(J,2));
%B=AdBd(1:size(B,1),end-size(B,2)+1:end);
Ad=(eye(N)+J*dt);
Bd=B*dt;

figure();
subplot(2,3,1);
latex_show('X',X);
subplot(2,3,4);
latex_show('U',U);
subplot(2,3,[2 3]);
latex_show('f',F);
subplot(2,3,[5 6]);
latex_show('h',h);

%% Check controllability and observability
C=Bd;
for i=1:N-1
    C=[C Ad^i*Bd];
end
if (rank(C)<N)
    warning(sprintf('System is not controllable rank(C)=%i<%i',rank(C),N))
end
O=H;
for i=1:N-1
    O=[O;H*Ad^i];
end
if (rank(O)<N)
    warning(sprintf('System is not observable rank(O)=%i<%i',rank(O),N))
end
%J=simplify(J);

%% Print as c code

warning off all
print_vec(F,'f',X,U);
print_mat(J,'J',X,U);
print_vec(h,'h',X,U);
print_mat(H,'H',X,U);
warning on all


