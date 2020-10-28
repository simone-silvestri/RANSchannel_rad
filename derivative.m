function dy = derivative(y,x)

n = length(x);
dy= zeros(n,1);

for i=2:n-1
    dx    = x(i+1)-x(i-1);
    dy(i) = (y(i+1)-y(i-1))/dx;
end
    dy(1)   = (y(2) - y(1))/(x(2)-x(1));
    dy(end) = (y(end) - y(end-1))/(x(end)-x(end-1));
end