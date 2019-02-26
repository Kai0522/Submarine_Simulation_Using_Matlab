function D=Damping_Matrix(x)%% Textbook Page 43
Xu=-9.30e-001;
Yv=-3.55e+001;
Ze=-3.55e+001;
Kp=-1.41e-002;
Mq=-4.88e+000;
Nr=-4.88e+000;

Xuu=-1.62e+000;
Yvv=-1.31e+002;
Zww=-1.31e+002;
Kpp=-1.30e-003;
Mqq=-9.40e+000;
Nrr=-9.40e+000;

D=-diag([Xu,Yv,Ze,Kp,Mq,Nr])-diag([Xuu*abs(x(1)),Yvv*abs(x(2)),Zww*abs(x(3)),Kpp*abs(x(4)),Mqq*abs(x(5)),Nrr*abs(x(6))]);
