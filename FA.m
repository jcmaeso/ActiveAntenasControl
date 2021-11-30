puntos = 1001;
u = linspace(-1,1,puntos);
v = linspace(-1,1,puntos);
savefilename = "prueba1";
nelm = 18;
separable = 1;
Ai = ones(1,nelm); %Amplitude distribution
Bi = ones(nelm,nelm); %Amplitude distribution non separable
alphay = 0;
alphax = 0;
sep = 0.5;
% H = 0.4; Pedestal
% n_inicio = -9.5;
% n_final = 9.5;
% n = linspace (n_inicio,n_final,nelm);
% for i = 1:nelm
%    Ai(i) = (1-H) + H*cos((pi*n(i))/(nelm-1))^2; Cosine in pedestal
% end
Ai = chebwin(nelm,20); %Chebwin
% Bi = Ai*Ai';


%Spatial Distribution
if mod(nelm,2) ~= 0
    xesp = ceil(-nelm/2)*sep:sep:floor(nelm/2)*sep;
else
    xesp = -nelm/2*sep+sep/2:sep:nelm/2*sep-sep/2;
end

[U,V] = meshgrid(u,v);

if separable == 1
    Fu = 0;
    for i = 1:nelm
        Fu = Fu + Ai(i).*exp(1j*(2*pi*xesp(i)*u+alphax*i)); 
    end
    Ai = taylorwin(nelm,2,-20);
    Fv = 0;
    for i = 1:nelm
        Fv = Fv + Ai(i).*exp(1j*(2*pi*xesp(i)*v+alphay*i)); 
    end
    F = Fu.'*Fv;
else
    F = 0;
    for i = 1:nelm
        for j = 1:nelm
            F = F + Bi(i,j).*exp(1j*((2*pi*xesp(i)*U+i*alphax)+(2*pi*xesp(j)*V+j*alphay))); 
        end
    end
end

phi = atand(U./V);
elementDiagram = zeros(puntos,puntos);

for i = 1:puntos
    for j = 1:puntos
        conv = hypot(u(i),v(j));
        if conv > 1
            continue;
        end
        elementDiagram(i,j) = cos(asin(conv)).^2;
    end 
end

diagram = elementDiagram.*F;
figure(1);
mesh(U,V,abs(diagram));
title("UV Plot")
save(savefilename,"u","v","F","elementDiagram","diagram");
%% Extract Phi Cut
[pat_azel,az_pat,el_pat] = uv2azelpat(abs(diagram),u,v,linspace(-90,90,puntos),linspace(-90,90,puntos));
figure(2);
plot(el_pat,20*log10(pat_azel(:,ceil(length(el_pat)/2)))-max(20*log10(pat_azel(:,ceil(length(el_pat)/2)))));
title("Phi = 0 Cut");
axis([-90,90,-40,0]);
figure(3);
plot(el_pat,20*log10(pat_azel(ceil(length(el_pat)/2),:))-max(20*log10(pat_azel(ceil(length(el_pat)/2),:))));
title("Theta = 0 Cut");
axis([-90,90,-40,0]);
