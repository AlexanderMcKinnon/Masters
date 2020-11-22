function dydt = Model(t, y, rates, parameters)
%% Unpacking
    %Parameters
    s0 = parameters(1);
	ns = parameters(2);
	nr = parameters(3);
	nt = parameters(4);
    nm = parameters(5);
    nq = parameters(6);
	gamma_max = parameters(7);
	Kp = parameters(8);
	vet = parameters(9);
	Ket = parameters(10);
	vem = parameters(11);
	Kem = parameters(12);
	wr_max = parameters(13);
	wt_max = parameters(14);
    wm_max = parameters(15);
	wq_max = parameters(16);
	thetar = parameters(17);
	thetat = parameters(18);
    thetam = parameters(19);
    thetaq = parameters(20);
	Kq = parameters(21);
    Nq = parameters(22);
	M = parameters(23);
    
    %Variables
    Mr = y(1);
	Mt = y(2);
	Mm = y(3);
	Mq = y(4);    
	Pr = y(5);
	Pt = y(6);
	Pm = y(7);
	Pq = y(8);        
    Cr = y(9);
	Ct = y(10);
	Cm = y(11);
	Cq = y(12);    
    Si = y(13);
	A = y(14);
    
    %Rates
    dm = rates(1);
	kb = rates(2);
	ku = rates(3);
    
    %Summations
    M_sum = Mr + Mt + Mm + Mq;
    C_sum = Cr + Ct + Cm + Cq;
%% Calculated Rates
    K_gamma = gamma_max/Kp;%unsure of this
    gamma = gamma_max*A/(K_gamma+A);
    ttrate = gamma*C_sum;
    lam = ttrate/M;
    vimp = Pt*vet*s0/(Ket+s0);
    vcat = Pm*vem*Si/(Kem+Si);
    
    vr = Cr*gamma/nr;
    vt = Ct*gamma/nt;
    vm = Cm*gamma/nm;
    vq = Cq*gamma/nq;
    
    wr = wr_max*A/(thetar + A);
    wt = wt_max*A/(thetat + A);
    wm = wm_max*A/(thetam + A);
    wq = (wq_max*A/(thetaq + A))*(1/(1+(Pq/Kq)^Nq));
%% ODES
    dydt(size(y,1),1) = 0;
    % mRNA (Ribo, Tran, Meta, HoKe)
    dydt(1) = wr - Mr*(lam+dm) - Pr*Mr*kb + Cr*ku + vr;
    dydt(2) = wt - Mt*(lam+dm) - Pr*Mt*kb + Ct*ku + vt;
    dydt(3) = wm - Mm*(lam+dm) - Pr*Mm*kb + Cm*ku + vm;
    dydt(4) = wq - Mq*(lam+dm) - Pr*Mq*kb + Cq*ku + vq;

    % Protein (Ribo, Tran, Meta, HoKe)
    dydt(5) = vr - Pr*lam + ku*C_sum - kb*Pr*M_sum + vr + vt + vq + vm;
    dydt(6) = vt - Pt*lam;
    dydt(7) = vm - Pm*lam;
    dydt(8) = vq - Pq*lam;

    % Complex (Ribo, Tran, Meta, HoKe)
    dydt(9) = kb*Pr*Mr - ku*Cr - vr - Cr*lam;
    dydt(10) = kb*Pr*Mt - ku*Ct - vt - Ct*lam;
    dydt(11) = kb*Pr*Mm - ku*Cm - vm - Cm*lam;
    dydt(12) = kb*Pr*Mq - ku*Cq - vq - Cq*lam;

    % Energy
    dydt(13) = vimp - vcat - Si*lam;
    dydt(14) = ns*vcat - (nr*vr) - (nt*vt) - (nm*vm) - (nq*vq) - lam*A;