%% Calculated Rates


%% ODES
% mRNA (Ribo, Tran, Meta, HoKe)
dydt(1) = wr - mr*(lam*dmr) - pr*mr*kb + Cr*ku + vr
dydt(2) = wt - mt*(lam*dmt) - pr*mt*kb + Ct*ku + vt
dydt(3) = wm - mm*(lam*dmm) - pr*mm*kb + Cm*ku + vm
dydt(4) = wq - mq*(lam*dmq) - pr*mq*kb + Cm*ku + vq

% Protein (Ribo, Tran, Meta, HoKe)
dydt(5) = vr - pr*lam + (- pr*mr*kb + Cr*ku + vr) + (- pr*mt*kb + Ct*ku + vt) + (- pr*mq*kb + Cm*ku + vq)
dydt(6) = vt - pt*(lam + dpt)
dydt(7) = vm - pm*(lam + dpm)
dydt(8) = vq - pq*(lam + dpq)

% Complex (Ribo, Tran, Meta, HoKe)
dydt(9) = kb*pr*mr - ku*Cr - vr - Cr*lam
dydt(10) = kb*pr*mt - ku*Ct - vt - Ct*lam
dydt(11) = kb*pr*mm - ku*Cm - vm - Cm*lam
dydt(12) = kb*pr*mq - ku*Cq - vq - Cq*lam

% Energy
dydt(13) = vimp - vcat - si*lam
dydt(14) = ns*vcat + (nr*vr) + (nt*vt) + (nm*vm) + (nq*vq) - lam*a