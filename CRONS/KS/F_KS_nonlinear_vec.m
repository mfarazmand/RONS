function Fij = F_KS_nonlinear_vec(x,q,r,param1)

q = reshape(q,[param1,r])';
q = repmat(q,[r,1]);

ql1=reshape(q(:,1),[r,r]);
ql2=reshape(q(:,2),[r,r]);
ql3=reshape(q(:,3),[r,r]);
ql4=reshape(q(:,4),[r,r]);

qj1=ql1';qj2=ql2';qj3=ql3';qj4=ql4';


innerpartj = qj2.*sin(x*pi/10 + qj3) + qj4;
innerpartl = ql2.*sin(x*pi/10 + ql3) + ql4;

Fij = -qj1.*ql1.*tanh(innerpartj).*...
    (ql2*(pi/10).*cos(x*pi/10 + ql3)).*sech(innerpartl).^2;

Fij = sum(Fij,'all');


end