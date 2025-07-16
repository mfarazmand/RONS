function F = fvec(q,x,r,param1,samples)

F = zeros(samples,1);
for ix = 1:samples
    F(ix) =  F_KS_linear_vec(x(ix),q,r,param1) +...
        F_KS_nonlinear_vec(x(ix),q,r,param1);
end