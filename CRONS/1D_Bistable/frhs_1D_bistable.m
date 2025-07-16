function dqdt = frhs_1D_bistable(t,q,r,param1,samples,delta,alpha1,alpha2,alpha3,nu)

t

dim = 1;
params = length(q);

x = linspace(-3,3,samples); % determine where to sample

M = zeros(samples,params); %metric tensor
F = zeros(samples,1); %rhs vector
gradI_m = zeros(params, 1);

for ix = 1:samples
    M(ix,:) = reshape(derivs_1D_bistable(x(ix),q,r,param1)',[1,params]);
    %easier to break  RHS into linear and nonlinear parts to vectorize each
    %of the terms
    F(ix) =  F0_1D_bistable_CRONS(x(ix),q,r,param1,alpha1,alpha2,alpha3,nu);
end


for i = 1:r
    gradI_m(param1*(i-1)+1:param1*i) = pi^(dim/2)*...
        [2*q(1+param1*(i-1)).*q(2+param1*(i-1)).^dim;
        dim*q(1+param1*(i-1)).^2.*q(2+param1*(i-1)).^(dim-1);
        zeros(dim,1)];
end


qk = reshape(q,[param1,r])';
q1 = qk(:,1); q2 = qk(:,2); q3 = qk(:,3);

denom = sum(q1.^2.*exp(-(x - q3).^2./q2.^2));
denom(denom<=1E-16)=1;

%% Do Fisher information metric
f = (M'./denom)*F;
M = (M'./denom)*M;

%% Calculate lambda to conserve total probability
minv= (M + delta * eye(params))\gradI_m;
finv= (M + delta * eye(params))\f;
lambda = dot(gradI_m, finv)/dot(gradI_m, minv);

%%
dqdt = (M + delta * eye(params)) \ (f- lambda*gradI_m);


end