%CRONS applied to Kuramoto--Sivashinsky equation
%To apply CRONS to another PDE with a different approximate solution, the
%files which calculate values of the metric tensor and RHS vectors need to be
%altered

clear
clc

%% input values

r = 10; % number of modes
param1 = 4; % number of parameters in one mode
samples = 128;
alpha = 1E-5; % regularization parameter

params = param1*r; % total number of parameters

%% LSQ fit for initial condition

% random starting values for parameters before fitting
A0 = -ones(1,r);
W0 = ones(1,r);
b0 = zeros(1,r);
XC0 = zeros(1,r);

x = linspace(-10,10,1000);

d = -sin(pi*x/10);

% q0 = param_estimation(q0,d,x,r); %least-squares fit

%fit for r = 10
q0 = [1.16816143930293,-0.0973029046515971,0.120602578325550,-0.000713715405178642,2.17057161427156,-0.109414890151566,-0.121022769214914,-1.39645166232610e-06,1.78442074698694,-0.154275629787411,0.0553338044612913,-1.02713730526624e-06,2.32223993159204,-0.0747700250963436,-0.0522880536682285,-1.11868729584162e-08,0.418058984849206,-0.157933363682466,0.184402070586593,-0.0686768615274412,0.288601653053627,0.402885520945006,0.0217978145987765,-0.00924024624338785,-0.215779248852785,0.179678315643780,0.0228225418794866,-0.179406615238103,1.37238494092780,-0.142751534558877,-0.0311237822644626,1.48400947087741e-07,0.00483493569461876,0.454340840069144,-0.544287789234260,-0.317717996473442,0.574032439346720,-0.0406890371082819,0.152777189109776,-0.00826090093871740];


%% time integration 

% time grid
t0 = 0.0;
tf = 40;
dt = 0.01;
t = t0 : dt : tf;

% rhs function
F = @(t, q)(frhs(t, q, r,samples,alpha,param1,params));

tic
options = odeset('RelTol',1E-4,'AbsTol',1E-5);
[t,Q] = ode45(F, t, q0,options);
% [t,Q] = ode113(F, t, q0,options);
% [t,Q] = ode15s(F, t, q0,options);
runtime = toc;

%% build solution
u = build_sol(x,t,Q,param1,r);

%% Plot solution
figure
pcolor(t,x,u)
shading interp
xlabel('$t$','interpreter','latex')
ylabel('$x$','interpreter','latex')
colormap('jet')
colorbar
shading interp
set(gca,'FontSize',20)


function dqdt = frhs(t,q,r,samples,alpha,param1,params)

t

x = linspace(-10,10,samples); % determine where to sample

M = zeros(samples,params); %metric tensor
f = zeros(samples,1); %rhs vector
for ix = 1:samples
    M(ix,:) = reshape(derivs(x(ix),q,r,param1)',[1,params]);
    %easier to break  RHS into linear and nonlinear parts to vectorize each
    %of the terms
    f(ix) =  F_KS_linear_vec(x(ix),q,r,param1) +...
        F_KS_nonlinear_vec(x(ix),q,r,param1);
end

dqdt = (M'*M + alpha * eye(params)) \ (M'*f); %regularized CRONS equation

end