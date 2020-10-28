function dy = derivativeX2(y,x)

n = length(x);
h = floor(n./2);
dy= zeros(n,1);
xu = zeros(n-1,1);
yu = zeros(n-1,1);

%creating xu
for i=1:n-1
    xu(i) = (x(i)+x(i+1))*0.5;
end

for i = 1:n-1
    yu(i) = (y(i+1)-y(i))./(x(i+1)-x(i));
end

dy(1:h) = interp1(xu(1:h-1),yu(1:h-1),x(1:h),'pchip');
dy(h+1:end) = interp1(xu(h:end),yu(h:end),x(h+1:end),'pchip');

end