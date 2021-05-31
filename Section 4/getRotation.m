function A = getRotation(d)
c = 1/sqrt(d)*ones(1,d);
A = [1/sqrt(d)*ones(1,d); null(c(:).')'];

