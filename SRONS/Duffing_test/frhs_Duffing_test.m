function dqdt = frhs_Duffing_test(t,q,r,param1,delta,alpha1,alpha2,alpha3,nu)

params = length(q);
M = zeros(params);
f = zeros(params, 1);

%build M using blocks
for i = 1:r
    for j = 1:i
        Mij = Mblock_Duffing_test(...
            [q(param1*(i-1)+1:param1*i),q(param1*(j-1)+1:param1*j)]);
        M(param1*(i-1)+1:param1*i,param1*(j-1)+1:param1*j) = Mij;
        M(param1*(j-1)+1:param1*(j),param1*(i-1)+1:param1*(i)) = Mij';
    end
end

qk = reshape(q,[param1,r])';

for k = 1:r
	f((1:param1) + param1*(k-1)) = ...
		F0block_Duffing_test(qk(k,:),q,r,param1,alpha1,alpha2,alpha3,nu);
end

dqdt = (M + delta * eye(params)) \ f;
 
end