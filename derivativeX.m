function dy = derivativeX(y,x)

n = length(x);
h = floor(n./2);
dy= zeros(n,1);
xu = zeros(n+1,1);
yu = zeros(n+1,1);

%creating xu
xu(1) = 0;
for i=2:n+1
    xu(i) = 2*x(i-1)-xu(i-1);
end

yu(1:h+1) = interp1(x(1:h),y(1:h),xu(1:h+1),'pchip');
yu(h+2:n+1) = interp1(x(h+1:n),y(h+1:n),xu(h+2:n+1),'pchip');

for i = 2:n+1
    dy(i-1) = (yu(i)-yu(i-1))./(xu(i)-xu(i-1));
end

end

