function Fij = F0_1D_bistable_CRONS(x,q,r,param1,alpha1,alpha2,alpha3,nu)

q  = reshape(q,[param1,r])';
q1 = q(:,1); q2 = q(:,2); q3 = q(:,3);

exppart = -(x - q3).^2./q2.^2;
uhat = q1.^2.*exp(exppart);
uhatx = (-2*(x - q3)./q2.^2).*uhat;
uhatxx = (-2./q2.^2).*uhat +...
    (-2*(x - q3)./q2.^2).*uhatx ;
Vx = 4*alpha1*x.^3 + 3*alpha2*x.^2 + 2*alpha3*x;
Vxx = 12*alpha1*x.^2 + 6*alpha2*x + 2*alpha3;

Fij = Vxx.*uhat + Vx.*uhatx + nu*uhatxx;

Fij = sum(Fij);
% Fij = sum(Fij./sum((uhat)));


end