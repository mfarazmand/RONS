function dqdt = frhs_duffing(t,q,r,param1,delta,alpha1,alpha2,alpha3,nu)

t

params = length(q);

M = zeros(params);
f = zeros(params, 1);
gradI_m = zeros(params, 1);  %used to enforce constant probability

qk = reshape(q,[param1,r])';


%% build M and f block by block
for i = 1:r
    for j = 1:i
        Mij = Mblock_Duffing([qk(i,:),qk(j,:)]);
        M(param1*(i-1)+1:param1*i,param1*(j-1)+1:param1*j) = Mij;
        M(param1*(j-1)+1:param1*(j),param1*(i-1)+1:param1*(i)) = Mij';
    end
    
    f((1:param1) + param1*(i-1)) = ...
        F0block_Duffing(qk(i,:),q,r,param1,alpha1,alpha2,alpha3,nu);
    
    gradI_m(param1*(i-1)+1:param1*i) = pi*...
    [2*q(1+param1*(i-1)).*q(2+param1*(i-1)).^2;
    2*q(1+param1*(i-1)).^2.*q(2+param1*(i-1));
    zeros(2,1)];
end

%% Calculate lambda value to enforce conserved quantity
minv= (M + delta * eye(params))\gradI_m;
finv= (M + delta * eye(params))\f;
lambda = dot(gradI_m, finv)/dot(gradI_m, minv);


%% Solve RONS equation
dqdt = (M + delta * eye(params)) \ (f- lambda*gradI_m);

end