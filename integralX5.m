function Iy = integralX5(y,x)

n = length(x);
Iy= zeros(n,1);
for i=2:n
    Iy(i) = Iy(i-1) + 0.5*(y(i)+y(i-1))*(x(i)-x(i-1));
end

end
