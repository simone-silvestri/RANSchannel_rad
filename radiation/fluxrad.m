function G = fluxrad(mesh,abs,E,Th,Tc)

abs = abs(1);
y = mesh.y;

n=length(y);

Ip = zeros(n,1);
Im = zeros(n,1);
Ip(1) = (Th./Tc).^4;
Im(n) = 1;

E = E./4;

for i=2:n
    dt = abs*(y(i)-y(i-1));
    Ib = 0.5*(E(i)+E(i-1));
    
    Ip(i) = 2*dt./(dt+1)*Ib - Ip(i-1)*(dt-1)/(dt+1);
end
for i=n-1:-1:1
    dt = -abs*(y(i)-y(i+1));
    Ib = 0.5*(E(i)+E(i+1));
    
    Im(i) = 2*dt./(dt+1)*Ib - Im(i+1)*(dt-1)/(dt+1);
end
for i=1:n
    G(i) = 2*(Ip(i)+Im(i));
end

end

