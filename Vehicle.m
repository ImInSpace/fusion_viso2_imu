%% State declaration
syms x          %x-coordinate
syms y          %y-coordinate
syms th         %yaw
syms th_dot     %yaw rate
syms Ve         %lateral velocity
syms Vn         %longitudinal velocity
X=[x;y;th;th_dot;Ve;Vn];

%% Input declaration
Pf = 0;          %Force to front tyres
syms Pr          %Force to rear tyres
syms d           %Steering angle
U=[Pr;d];

%% Parameter declaration
%m=1292.2;       %Vehicle mass
%I=2380.7;       %Vehicle inertia
%a=1.006;        %CG to front axle
%b=1.534;        %CG to rear axle
syms m I a b
h=0.3;          %CG to ground
g=9.81;         %gravity     
mu=0.85;        %friction coeff
c=20000;        %cornering stiffness

%% Segel Lateral Force Model
ignore_lateral_force=1;
if ~ignore_lateral_force
    Fzf=(m*g*b-(Pf+Pr)*h)/(a+b); %normal load on front tyre
    Fzr=(m*g*a-(Pf+Pr)*h)/(a+b); %normal load on rear tyre

    af=d-(a*th_dot+Ve)/Vn;  %front wheel slip angle
    ar=(b*th_dot-Ve)/Vn;    %rear wheel slip angle

    asf=c*af/(mu*Fzf);
    asr=c*ar/(mu*Fzr); 

    Fef=mu*Fzf*(asf-(asf*abs(asf))/3+asf^3/27)*sqrt(1-Pf^2/(mu^2*Fzf^2)+Pf^2/c^2); %Lateral force on front tyre
    Fer=mu*Fzr*(asr-(asr*abs(asr))/3+asr^3/27)*sqrt(1-Pr^2/(mu^2*Fzr^2)+Pr^2/c^2); %Lateral force on rear tyre
else
    Fef=Pr*d/(1+d);
    Fer=-Pr*d/(1+d);
end

%% ODE definition
F=[-Ve*sin(th)+Vn*cos(th);
    Ve*cos(th)+Vn*sin(th);
    th_dot;
    (a*Pf*d+b*Fef-b*Fer)/I;
    (Pf*d+Fef+Fer)/m-Vn*th_dot;
    (Pf+Pr-Fef*d)/m-Vn*th_dot;
  ];

F=simplify(F);
J=jacobian(F,X);
J=simplify(J);


syms x(t) u(t)
F=subs(F,X,x(0:length(X)-1).');
F=subs(F,U,u(0:length(U)-1).');
J=subs(J,X,x(0:length(X)-1).');
J=subs(J,U,u(0:length(U)-1).');
fprintf("f << ");
for i=1:6
    Fi=ccode(F(i));
    Fi=regexprep(Fi,'  t0 = ','');
    Fi=regexprep(Fi,';','');
    Fi=regexprep(Fi,'.0)',')');
    if(i==6)
        fprintf('%s;\n\n',Fi)
    else
        fprintf('%s,\n',Fi)
    end
end

for i=1:6
    for j=1:6
        if J(i,j)~=0
            Jij=ccode(J(i,j));
            Jij=regexprep(Jij,'  t0 = ','');
            Jij=regexprep(Jij,';','');
            Jij=regexprep(Jij,'.0)',')');
            fprintf('J(%i,%i)=%s;\n',i-1,j-1,Jij)
        end
    end
end