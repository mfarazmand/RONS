function dudq = derivs(x,q,r,param1)

%reshape parameters into matrix where first column is amplitudes,
%second column is beta1,  3rd column beta2, etc
q = reshape(q,[param1,r])';


% define derivs with respect to each parameter
innerpart = q(:,2).*sin(x*pi/10+q(:,3)) + q(:,4);
dudq = [tanh(innerpart),...
    q(:,1).*sin(x*pi/10+q(:,3)).*sech(innerpart).^2,...
    q(:,1).*q(:,2).*cos(x*pi/10+q(:,3)).*sech(innerpart).^2,...
    q(:,1).*sech(innerpart).^2];

end