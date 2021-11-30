clc;
clear;

nelm = 20;
sep = 210;
pos =  [];
if mod(nelm,2) ~= 0
    pos = ceil(-nelm/2)*sep:sep:floor(nelm/2)*sep;
else
    pos = -nelm/2*sep+sep/2:sep:nelm/2*sep-sep/2;
end


phax = 0*(1:nelm);
phay = 0*(1:nelm);
%Ai1 = taylorwin(nelm,3,-25); %Chebwin
Ai1 = ones(nelm,1);
Ai2 = taylorwin(nelm,3,-20); %Chebwin
H = 0.4;
n_inicio = -9.5;
n_final = 9.5;
n = linspace (n_inicio,n_final,nelm);
for i = 1:nelm
   Ai2(i) = (1-H) + H*cos((pi*n(i))/(nelm-1))^2; %Cosine in pedestal
end
Bi = Ai1*Ai2';
Bi = ones(nelm,nelm);

cst = actxserver('CSTStudio.Application');
mws = invoke(cst, 'Active3D');
%SelectTreeItem = mws.invoke('SelectTreeItem',"Farfields\farfield (f=60) [1]");
% farfield_plot = invoke(mws,'FarfieldPlot');
% invoke(farfield_plot,'Plottype','cartesian');
% invoke(farfield_plot,'Step',0.1);
% invoke(farfield_plot,'SetAntennaType','directional_linear');
% invoke(farfield_plot,'SetTheta360',true);
% invoke(farfield_plot,'SetPlotRangeOnly',true);
% invoke(farfield_plot,'Plot');

SelectTreeItem = mws.invoke('SelectTreeItem',"Farfields\farfield (f=1.02) [1(2)]");
farfield_ctr = invoke(mws,'FarfieldArray');
invoke(farfield_ctr,'Reset');
invoke(farfield_ctr,'UseArray',true);
invoke(farfield_ctr,'Arraytype',"edit");
% for i = 1:nelm
%     for j = 1:nelm
%         invoke(farfield_ctr,'Antenna',pos(i),pos(j),0,Bi(i,j),phax(i)+phay(j));
%     end
% end
for i = 1:nelm
    %Los dos primeros te dan igual
    %
  invoke(farfield_ctr,'Antenna',pos(i),0,0,Ai2(i),phay(i));
end

invoke(farfield_ctr,'SetList');

%SelectTreeItem = mws.invoke('SelectTreeItem',"Farfields\farfield (f=60) [1]\Ludwig 3 Crosspolar");







