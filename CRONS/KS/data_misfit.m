function J = data_misfit(p,d,x,r)
%   p --- parameter vector
%   d --- vector of measurement data
y = zeros(1,length(x));
for iq = 0:r-1
    y = y + p(1+iq*param1)*...
        tanh(p(2+iq*param1)*sin(x*pi/10 + p(3+iq*param1)) + p(4+iq*param1));
end
% for iq = 0:r-1
%     y = y + p(1+iq*4)*tanh(p(2+iq*4)*sin(x*pi/10 + p(3+iq*4)) + p(4+iq*4));
% end
% J = norm(z-d).^2;

J = norm(y-d)

