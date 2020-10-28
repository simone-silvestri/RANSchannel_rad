function I = integralX(x,y)

I    = zeros(length(y),1);

for i=2:length(x)
    I(i) = I(i-1) + 0.5*(y(i)+y(i-1))*(x(i)-x(i-1));
end

end

