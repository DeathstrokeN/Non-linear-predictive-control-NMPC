function J = evalcrit(Hp, Hc, T, lambda, x, rrad, u, ukm1, K, m, L, g) 
rpred = rrad*ones(Hp,1);
%%
%vector deltau
deltau = zeros(Hc,1);
deltau(1) = u(1) - ukm1;
if Hc>1   
    for i=2:Hc
        deltau(i) = u(i) - u(i-1);
    end
end
%%
Ba = [0; T/(m*L^2); 0];
%predictions ypred
%ypred = x(1)+x(3)...(xe(1)+xi)
ypred = zeros(Hp,1);
for i=1:Hc
    ypred(i) = T*x(2)+x(1)+x(3);
    x = [(T*x(2)+x(1)); ((-g*T/L)*sin(x(1))-(K*T/m)*x(2)+x(2)); x(3)] + Ba*u(i);    
end
for i=Hc+1:Hp
    ypred(i) = T*x(2)+x(1)+x(3);
    x = [(T*x(2)+x(1)); ((-g*T/L)*sin(x(1))-(K*T/m)*x(2)+x(2)); x(3)]+ Ba*u(Hc);    
end

%%
%criterion
J = sum((ypred - rpred).^2) + lambda*sum(deltau.^2);
end
