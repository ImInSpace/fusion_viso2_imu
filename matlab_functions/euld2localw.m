function euld2localw=euld2localw(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    euld2localw=[-sth 0 1;
    cth*sph cph 0;
    cph*cth -sph 0];
end