clear all
close all

s0 = 10000;
ns = 0.5;
nr = 7549;
nt = 300;
nm = 300;
nq = 300;
gamma_max = 1260;
Kp = 180.1;
vet = 726;
Ket = 1000;
vem = 5800;
Kem = 1000;
wr_max = 930;
wt_max = 4.139;
wm_max = 4.139;
wq_max = 948.9;
thetar = 426.9;
thetat = 4.38;
thetam = 4.38;
thetaq = 4.38;
Kq = 152220;
Nq = 4;
M = 1e8;
parameters = [s0, ns, nr, nt, nm, nq, gamma_max, Kp, vet, Ket, vem, Kem, wr_max, wt_max, wm_max, wq_max, thetar, thetat, thetam, thetaq, Kq,Nq, M];

pr_0 = 10;
pt_0 = 0;
pm_0 = 0;
pq_0 = 0;    
mr_0 = 0;
mt_0 = 0;
mm_0 = 0;
mq_0 = 0;        
cr_0 = 0;
ct_0 = 0;
cm_0 = 0;
cq_0 = 0;    
si_0 = 0;
a_0 = 1000;
init = [mr_0, mt_0, mm_0, mq_0, pr_0, pt_0, pm_0, pq_0, cr_0, ct_0, cm_0, cq_0, si_0, a_0];

dm = 0.1;
kb = 1;
ku = 1;
rates = [dm, kb, ku];

t0 = 0;
tf = 600;
[t,y] = ode15s(@(t,y) Model(t, y, rates, parameters), [t0 tf], init);

Mr = y(:,1);
Mt = y(:,2);
Mm = y(:,3);
Mq = y(:,4);
Pr = y(:,5);
Pt = y(:,6);
Pm = y(:,7);
Pq = y(:,8);
Cr = y(:,9);
Ct = y(:,10);
Cm = y(:,11);
Cq = y(:,12);
Si = y(:,13);
A = y(:,14);

figure(1)
subplot(2,2,1)
hold on
title("Proteins"); ylabel("Magnitude"); xlabel("time")
plot(t,Pr,'r')
plot(t,Pt,'g')
plot(t,Pm,'b--')
plot(t,Pq,'c')
legend(["Ribo","Tran","Meta","HoKe"])

subplot(2,2,2)
hold on
title("mRNA"); ylabel("Magnitude"); xlabel("time")
plot(t,Mr,'r')
plot(t,Mt,'g')
plot(t,Mm,'b--')
plot(t,Mq,'c')
legend(["Ribo","Tran","Meta","HoKe"])

subplot(2,2,3)
hold on
title("Complexes"); ylabel("Magnitude"); xlabel("time")
plot(t,Cr,'r')
plot(t,Ct,'g')
plot(t,Cm,'b--')
plot(t,Cq,'c')
legend(["Ribo","Tran","Meta","HoKe"])

subplot(2,2,4)
hold on
title("Energy"); ylabel("Magnitude"); xlabel("time")
plot(t,Si,'y')
plot(t,A,'m')
legend(["InNut","Ener"])

figure(2)
subplot(2,2,1)
hold on
title("Ribosomes"); ylabel("Magnitude"); xlabel("time")
plot(t,Pr,'r')
plot(t,Mr,'g')
plot(t,Cr,'b')
legend(["Pro","mRNA","Comp"])

subplot(2,2,2)
hold on
title("Transport Proteins"); ylabel("Magnitude"); xlabel("time")
plot(t,Pt,'r')
plot(t,Mt,'g')
plot(t,Ct,'b')
legend(["Pro","mRNA","Comp"])

subplot(2,2,3)
hold on
title("Metabolic Proteins"); ylabel("Magnitude"); xlabel("time")
plot(t,Pm,'r')
plot(t,Mm,'g')
plot(t,Cm,'b')
legend(["Pro","mRNA","Comp"])

subplot(2,2,4)
hold on
title("House-Keeping Proteins"); ylabel("Magnitude"); xlabel("time")
plot(t,Pq,'r')
plot(t,Mq,'g')
plot(t,Cq,'b')
legend(["Pro","mRNA","Comp"])