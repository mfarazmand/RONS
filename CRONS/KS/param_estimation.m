function pstar = param_estimation(p0,d,x,r)
% purpose: parameter estimation routine
%
% input: 
%    p0  --- initial iterate to be used by optimization routine
%    d   --- vector of measurement data
%
% NOTE: this function requires the function "data_misfit"
%       that needs to be written as a part of the homework 
%       assignment. 

% estimate parameters by minimizing the data misfit function
J     = @(p)(data_misfit(p,d,x,r));
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
options = optimoptions('fmincon','MaxFunctionEvaluations',1E5);
pstar = fmincon(J,p0,A,b,Aeq,beq,lb,ub,[],options);
% pstar = fmincon(J,p0);



% fprintf('estimated parameters: (A, L, V,phi) = ...(%6.4f, %6.4f, %6.4f %6.4f)\n', ...
%     A, L, V, phi);
