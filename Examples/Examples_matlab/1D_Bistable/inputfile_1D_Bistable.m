clc
clear

r = 10; %number of modes
dim = 1; %dimension of problem
param1 = 2 + dim; %number of params in one mode
delta = 1E-4; %regularization parameter

A0 = ones(1,r)/sqrt(20);
L0 = 1./( sqrt(pi)/2 )*ones(1,r);
XC0 = [-ones(1,r/2), -2*ones(1,r/2)];

% A0 = rand(1,r);
% L0 =  rand(1,r);
% XC0 = 4*(rand(1,r));

volans = sqrt(pi)*sum(A0.^2.*L0); %used to normalize

q0=zeros(r*param1,1); %Initial condition
for i = 0:r-1
    q0((1:param1)+i*param1) = ...
        [A0(i+1)/sqrt(volans),L0(i+1),XC0(i+1)];
end

sigma = 0.5; %noise level
nu = sigma^2/2;

alpha1 = 1/4;
alpha2 = 0;
% alpha2 = -1/27;
alpha3 = -1/2;


t0 = 0.0;
ts = 5; %when to speed up
tf = 400;
% t = t0 : 1 : tf;
dt1 = 0.1;
dt2 = 2.5;
t = [t0:dt1:ts - dt1,...
    ts:dt2:tf];


% rhs function
F = @(t, q)(frhs_1D_bistable(t, q,r,param1,delta,alpha1,alpha2,alpha3,nu));

% Call solver
tic
options = odeset('RelTol',1e-6,'AbsTol',1e-6); 
% [t,Q] = ode45(F, t, q0,options);
% [t,Q] = ode113(F, t, q0,options);
[t,Q] = ode15s(F, t, q0,options);
runtime = toc;

nt = length(t);


x = -3:.01:3;
u = zeros(length(t),length(x));

for j = 1:nt
    for iq = 1:r
        gaussnum = (iq-1)*param1;
        u(j,:,:)  = squeeze(u(j,:,:)) + Q(j,1+gaussnum).^2.*...
            exp(-( (x - Q(j,3+gaussnum)).^2  )./Q(j,2+gaussnum).^2);
    end
     disp(['t=',num2str(t(j))]);
end



%% Analytical equilibrium
B = @(x) alpha1*x.^4 + alpha2*x.^3 + alpha3*x.^2;
Bprime = @(x) 4*alpha1*x.^3 + 3*alpha2*x.^2 + 2*alpha3*x;


vol = integral(@(x)exp( -B(x)/nu ),...
        -Inf,Inf,'RelTol',0,'AbsTol',1e-12);

uF = exp( -B(x)/nu )/vol;


%% Calculate total probability
mass = nan(1,length(t));
for j = 1:length(t)
    mass(j)  = sqrt(pi)*sum(Q(j,1:param1:end).^2.*Q(j,2:param1:end));
end

%% Equilibrium error
figure
plot(x,u(end,:) - uF,'linewidth',1)
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$\hat{p}_{eq} - p_{eq}$','interpreter','latex','fontsize',14)
title('Equilibrium Error' ,'interpreter','latex','fontsize',14)
set(gca,'FontSize',20)

%% Movie
M = struct('cdata', cell(1,ceil(length(t))), 'colormap', cell(1,ceil(length(t)))); %intialize movie
mcount = 1;
fig2 =  figure;
for it = 1:length(t)
    box on
    plot(x, uF,'k--','linewidth',2)
    hold on
    plot(x, u(it,:),'b','linewidth',2)
    hold off
    xlabel('$x$','interpreter','latex')
    ylabel('$p$','interpreter','latex')
    title(['PDF, $t = ',num2str(t(it),'%4.2f'),'$' ],'interpreter','latex','fontsize',14)
    ylim([0, max(max(u))+0.01])
    xlim([x(1), x(end)])
%     axis tight
    set(gca,'FontSize',24)
         set(gcf,'color','w');
    M(mcount) = getframe(fig2);
    mcount = mcount+1;
end

