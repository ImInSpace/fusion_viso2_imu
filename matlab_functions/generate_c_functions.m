
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

do_plot=0
if do_plot
    figure();
    subplot(2,3,1);
    latex_show('\v x',X);
    subplot(2,3,4);
    latex_show('\v u',U);
    subplot(2,3,[2 3]);
    latex_show('f(\v x,\v u)',F);
    subplot(2,3,[5 6]);
    latex_show('h(\v x)',h);
else
    sep=' \quad ';
    clipboard('copy',[latex_show('\v x',X),sep,...
    latex_show('\v u',U),sep,...
    latex_show('f(\v x,\v u)',F),sep,...
    latex_show('h(\v x)',h)...
    ])
end

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
print_vec(F,'f',X,U,1);
print_mat(J,'J',X,U,1);
print_vec(h,'h',X,U,1);
print_mat(H,'H',X,U,1);
warning on all


