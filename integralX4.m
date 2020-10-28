function Iy = integralX4(y,x)

n = length(x);
Iy= zeros(n,1);
yw = interp1(x,y,0.0,'pchip');
Iy(1) = 0.5*(y(1)+yw)*x(1);
for i=2:n
    Iy(i) = Iy(i-1) + 0.5*(y(i)+y(i-1))*(x(i)-x(i-1));
end

end
