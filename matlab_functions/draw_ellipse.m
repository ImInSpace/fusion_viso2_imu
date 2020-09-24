function ellipse = draw_ellipse(X,Q)
    [V,D]=eig(Q);
    radii=sqrt(diag(D)*chi2inv(0.95,2));
    th=linspace(0,2*pi,20);
    th(end)=[];
    x=radii(1)*sin(th);
    y=radii(2)*cos(th);
    P=X(:)+V*[x;y];
    ellipse = polyshape(P(1,:),P(2,:));
end