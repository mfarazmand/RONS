function u = build_sol(x,t,Q,param1,r)

u = zeros(length(x), length(t));
count = 1;
for it = 1:length(t)
    for iq = 0:r-1
        u(:,count) = u(:,count)' + Q(it,1+iq*param1)*...
            tanh(Q(it,2+iq*param1)*sin(x*pi/10 + Q(it,3+iq*param1) )+Q(it,4+iq*param1));
    end
    count = count+1;
end