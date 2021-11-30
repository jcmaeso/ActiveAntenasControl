nelm = 8;
separation = 0.55;
test_antena = ActiveAntenna(nelm,separation,6,0.5,true);

%% Phases Test
aims = [5,10,15,30];
phases_real = zeros(length(aims),nelm);
phases_dis = zeros(length(aims),nelm);
for i = 1:length(aims)
    [phases_dis(i,:),phases_real(i,:),~] = test_antena.aim2dir(aims(i));
end

figure(1);
for i = 1:length(aims)
    [diagram,theta] = array_factor(nelm,separation,phases_dis(i,:),ones(1,nelm));
    plot(theta,diagram);
    if(i == 1)
        hold on;
    end
end
hold off;
xlim([-45,45]);
ylim([-40,0]);
legend(["5","10","15","30"]);

%% Amplitudes Test
amps_real = zeros(5,nelm);
amps_dis = zeros(5,nelm);
titles = ["Coseno sobre pedestal","Cuadrada","Triangular","Chebyshev","Taylor"];
[amps_dis(1,:),amps_real(1,:)] = test_antena.amps_cosine_over_pedestal(0.4,10);
[amps_dis(2,:),amps_real(2,:)] = test_antena.amps_parabolic(0.5);
[amps_dis(3,:),amps_real(3,:)] = test_antena.amps_triangular(2);
[amps_dis(4,:),amps_real(4,:)] = test_antena.amps_cheb(20);
[amps_dis(5,:),amps_real(5,:)] = test_antena.amps_taylor(20,5);

figure(2);
for i = 1:5
    subplot(3,2,i)
    plot(1:nelm,amps_real(i,:),"*");
    hold on;
    plot(1:nelm,amps_dis(i,:),"*");
    title(titles(i));
    hold off;
    xlim([1,nelm]);
end

figure(3);
for i = 1:5
    subplot(3,2,i)
    [diagram,theta] = array_factor(nelm,separation,zeros(1,nelm),10.^(amps_real(i,:)/20));
    plot(theta,diagram);
    hold on;
    [diagram,theta] = array_factor(nelm,separation,zeros(1,nelm),10.^(amps_dis(i,:)/20));
    plot(theta,diagram);
    xlim([-45,45]);
    ylim([-40,0]);
    title(titles(i));
    hold off;
end

