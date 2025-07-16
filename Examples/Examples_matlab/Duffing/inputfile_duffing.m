clc
clear

r = 30; %number of modes
dim = 2; %dimension of problem
param1 = 2 + dim; %number of params in a mode
delta = 1E-3; %regularization parameter

A0 = 1*ones(1,r);
L0 = 1./( sqrt(r*pi*A0.^2) );
XC0 = [-ones(1,r/2), ones(1,r/2)];
YC0 = [-ones(1,r/2), ones(1,r/2)];

volans = pi*sum(A0.^2.*L0.^2); %used to normalize

for i = 0:r-1
    q0([1+i*param1,2+i*param1,3+i*param1,4+i*param1]) = ...
        [A0(i+1)/sqrt(volans),L0(i+1),XC0(i+1),YC0(i+1)];
end

%% parameter values
sigma = sqrt(1/20); %noise level
nu = sigma^2/2;

%PDE params
alpha1 = 1;
alpha2 = -0.2;
alpha3 = -1;

%% Time integration
% time grid
t0 = 0.0;
tf =  50;
% t = t0 : 0.5 : tf;
t = [0 : 0.1 : 25,25.5:0.5:50];


% rhs function
F = @(t, q)(frhs_duffing(t, q,r,param1,delta,alpha1,alpha2,alpha3,nu));

% Call solver
tic
options = odeset('RelTol',1e-3,'AbsTol',1e-6); 
% options = odeset('RelTol',1e-8,'AbsTol',1e-4); 
[t,Q] = ode45(F, t, q0,options);
% [t,Q] = ode113(F, t, q0,options);
% [t,Q] = ode15s(F, t, q0,options);
runtime = toc;
nt = length(t);

%% Build solution

x = -2:.01:2;
y = x';
u = zeros(length(t), length(y),length(x));

for j = 1:nt
    for iq = 1:r
        gaussnum = (iq-1)*param1;
        u(j,:,:)  = squeeze(u(j,:,:)) + Q(j,1+gaussnum).^2.*exp(-( (x - Q(j,3+gaussnum)).^2 +...
            (y - Q(j,4+gaussnum)).^2 )./Q(j,2+gaussnum).^2);
    end
     disp(['t=',num2str(t(j))]);
end

%% Analytical equilibrium
V = @(x,y) 4*x.^2 - 2*x.^4 - 4*y.^2;
vol = integral2(@(x,y)exp( V(x,y) ),...
        -Inf,Inf,-Inf,Inf,'RelTol',0,'AbsTol',1e-12);
uF = exp( V(x,y) )/vol;

%% Calculate total probability
mass = nan(1,length(t));

for j = 1:length(t)
    mass(j)  = pi*sum(Q(j,1:param1:end).^2.*Q(j,2:param1:end).^2);
end

%% Plot Equilibriums and error

figure
pcolor(x,y,uF)
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y$','interpreter','latex','fontsize',14)
zlabel('$p$','interpreter','latex','fontsize',14)
xlim([-2,2])
ylim([-2,2])
title('$p_{eq}$','interpreter','latex','fontsize',14)
colorbar
% caxis([cmin cmax])
shading interp
set(gca,'FontSize',24)

figure
pcolor(x,y,squeeze(u(end,:,:)) )
shading interp
colorbar
xlim([-2, 2])
ylim([-2, 2])
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y$','interpreter','latex','fontsize',14)
xlim([-2,2])
ylim([-2,2])
title('$\hat{p}_{eq}$','interpreter','latex','fontsize',14)
colorbar
% caxis([cmin cmax])
shading interp
set(gca,'FontSize',24)

figure
pcolor(x,y,squeeze(u(end,:,:)) - uF)
shading interp
colorbar
xlim([-2, 2])
ylim([-2, 2])
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y$','interpreter','latex','fontsize',14)
title('$\hat{p}_{eq} - p_{eq}$','interpreter','latex','fontsize',14)
set(gca,'FontSize',24)


%% Movie

Mmovie = struct('cdata', cell(1,ceil(length(t))), 'colormap', cell(1,ceil(length(t)))); %intialize movie
mcount = 1;
fig2 =  figure;
for it = 1:1:length(t)
    pcolor(x,y,squeeze(u(it,:,:)))
    shading interp
    hold off
    xlabel('$x$','interpreter','latex')
    ylabel('$y$','interpreter','latex')
    title(['PDF, $t = ',num2str(t(it),'%4.2f'),'$' ],'interpreter','latex','fontsize',14)
    ylim([-2., 2.])
    xlim([-2., 2.])
    zlim([0, 1])
%     ylim([0, 0.9])
%     caxis([0, 0.8])
%     axis tight
    axis tight
    axis square
    colormap(jet)
    colorbar
    set(gca,'FontSize',24)
         set(gcf,'color','w');
    Mmovie(mcount) = getframe(fig2);
    mcount = mcount+1;
end
