function lw2euld=localw2euld(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    lw2euld=[0 sph/cth cph/cth;
    0 cph -sph;
    1 sph*sth/cth cph*sth/cth];
end