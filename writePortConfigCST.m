function [diagrams] = writePortConfigCST(file,amps,phases,labels)
%WRITEPORTCONFIGCST Summary of this function goes here
%   Detailed explanation goes here
n_configs = size(amps,1); %Get number of iterations
%Open CST file
cst = actxserver('CSTStudio.Application');
mws = invoke(cst, 'OpenFile', fullfile(cd,file));
for i = 1:n_configs
    cmb = mws.invoke('CombineResults');
    cmb.invoke('Reset');
    cmb.invoke('SetMonitorType','frequency')
    cmb.invoke('EnableAutomaticLabeling','False')
    cmb.invoke('SetLabel',labels(i))
    cmb.invoke('SetNone')
    for j = 1:size(amps,2)
        cmb.invoke("SetExcitationValues","port",sprintf("%d",j),1,amps(i,j),phases(i,j));
    end    
    cmb.invoke('Run')
    
    %Default mode 1 no multimode taken into account at port
end
close(cmb);
mws.invoke("save");
mws.invoke("quit");
end
