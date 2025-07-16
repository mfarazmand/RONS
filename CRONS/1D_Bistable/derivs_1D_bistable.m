function dudq = derivs_1D_bistable(x,q,r,param1)

%reshape parameters into matrix where first column is amplitudes,
%second column is beta1,  3rd column beta2, etc
q = reshape(q,[param1,r])';
q1 = q(:,1); q2 = q(:,2); q3 = q(:,3);


% define derivs with respect to each parameter
exppart = exp(-(x - q3).^2./q2.^2);
uhat = q1.^2.*exppart;
dudq = [2*q1.*exppart,...
    2*(x - q3).^2./q2.^3.*uhat,...
    2*(x - q3)./q2.^2.*uhat];
% dudq = [2*q1.*exppart./sum(sqrt(uhat)),...
%     2*(x - q3).^2./q2.^3.*uhat./sum(sqrt(uhat)),...
%     2*(x - q3)./q2.^2.*uhat./sum(sqrt(uhat))];

end