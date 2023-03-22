% P is the product. N is the substrate, i denotes the rare isotope
step = 3e5;
%% reservoir 1
P = zeros(1,step);
Pi = zeros(1,step);
dP = zeros(1,step);
dPi = zeros(1,step);
N = zeros(1,step);
Ni = zeros(1,step);

Pdelta = zeros(1,step);
dPdelta = zeros(1,step);
Ndelta = zeros(1,step);
R_ref = 0.004;
delta_15 = 50;
Ndelta(1) = delta_15;
N(1) = 1e2;
Ni(1) = ( delta_15/1000 + 1 )*R_ref*N(1);
%% reservoir 2
P2 = zeros(1,step);
Pi2 = zeros(1,step);
dP2 = zeros(1,step);
dPi2 = zeros(1,step);
N2 = zeros(1,step);
Ni2 = zeros(1,2.5e4);

Pdelta2 = zeros(1,step);
dPdelta2 = zeros(1,step);
Ndelta2 = zeros(1,step);
delta_15_2 = 5;
Ndelta2(1) = delta_15_2;
N2(1) = 2e1;
Ni2(1) = ( delta_15_2/1000 + 1 )*R_ref*N2(1);

e = 10;
alpha = 1 - e/1000;
k = 0.2e-2;
ki = k*alpha;

dt = 0.01;

%%
N_in = 1e-1;
delta_in = 10;
Ni_in = ( delta_in./1e3 + 1 ).*R_ref.*N_in;
%%
for I = 2:step
    P(I) = P(I-1) + k.*N(I-1).*dt;
    dP(I) = k.*N(I-1).*dt;
    N(I) = N(I-1) - ( P(I) - P(I-1) ) + N_in.*dt;  
    Pi(I) = Pi(I-1) + ki.*Ni(I-1).*dt;
    dPi(I) = ki.*Ni(I-1).*dt;
    Ni(I) = Ni(I-1) - ( Pi(I) - Pi(I-1) ) + Ni_in.*dt;
    Pdelta(I) = ( Pi(I)./P(I)./R_ref - 1).*1e3;
    dPdelta(I) = ( dPi(I)./dP(I)./R_ref - 1).*1e3;
    Ndelta(I) = ( Ni(I)./N(I)./R_ref - 1).*1e3;
    
    P2(I) = P2(I-1) + k.*N2(I-1).*dt;
    dP2(I) = k.*N2(I-1).*dt;
    N2(I) = N2(I-1) - ( P2(I) - P2(I-1) ) + N_in.*dt;  
    Pi2(I) = Pi2(I-1) + ki.*Ni2(I-1).*dt;
    dPi2(I) = ki.*Ni2(I-1).*dt;
    Ni2(I) = Ni2(I-1) - ( Pi2(I) - Pi2(I-1) ) + Ni_in.*dt;
    Pdelta2(I) = ( Pi2(I)./P2(I)./R_ref - 1).*1e3;
    dPdelta2(I) = ( dPi2(I)./dP2(I)./R_ref - 1).*1e3;
    Ndelta2(I) = ( Ni2(I)./N2(I)./R_ref - 1).*1e3;
    
    
end

figure(1);
plot(N,"LineWidth",2); hold on;
plot(N2,"LineWidth",2); hold on;
xlimit = get(gca,'Xlim');
line(xlimit,[N(end) N(end)],'Color','black','LineStyle',':','LineWidth',1); hold on;
txt1 = ['k =' num2str(k)];
t1 = text(2e5,30,txt1);
t1.FontSize = 17;
txt2 = ['N_{input} =' num2str(N_in)];
t2 = text(2e5,40,txt2);
t2.FontSize = 17;
xlabel('time');
ylabel('substrate amount');
legend("reservoir 1: N_{0} = 100","reservoir 2: N_{0} = 20");
set(gca,'FontSize',14,'linewidth',1);
print(gcf,'-r600','-djpeg',['t_N','.jpeg']);

figure(2);
plot(Ndelta,"LineWidth",2); hold on;
plot(Ndelta2,"LineWidth",2); hold on;
xlimit = get(gca,'Xlim');
line(xlimit,[Ndelta(end) Ndelta(end)],'Color','black','LineStyle',':','LineWidth',1); hold on;
txt = ['\epsilon =' num2str(e)];
t = text(2e5,10,txt);
t.FontSize = 17;
txt3 = ['\delta^{15}N_{input} =' num2str(delta_in) '‰'];
t3 = text(2e5,15,txt3);
t3.FontSize = 17;
xlabel('time');
ylabel('substrate \delta^{15}N');
legend("reservoir 1: \delta^{15}N_{0} = 50‰","reservoir 2: \delta^{15}N_{0} = 5‰");
set(gca,'FontSize',14,'linewidth',1);
print(gcf,'-r600','-djpeg',['t_delta','.jpeg']);


    