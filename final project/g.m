% % simple person to person model
% % % Q=fg(r,r0,r1)t;Q=1
% % % f is a constant[1,400](for hour), g(r) is function of distance r, t is time(hour)
% % % g(r,r0,r1)=1/r2 when r0<=r<=r1 , g(r)=1 when r<r0, g(r)=0 when r>r1
% % % 

function Q=g(r,r0,r1)
if r<r0
    Q=1;
elseif r<=r1
   Q=1/r;
else
   Q=0;
end
end
