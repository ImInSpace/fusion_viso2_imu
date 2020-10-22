function rotm=my_eul2rotm(eul)
    psi=eul(1);theta=eul(2);phi=eul(3);
    cps=cos(psi);sps=sin(psi);
    cth=cos(theta);sth=sin(theta);
    cph=cos(phi);sph=sin(phi);
    rotm=[cps*cth cps*sph*sth-cph*sps sph*sps+cph*cps*sth;
        cth*sps cph*cps+sph*sps*sth cph*sps*sth-cps*sph;
        -sth cth*sph cph*cth];
end