%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read me %%%%%%%%%%%%%%%%%%%%%%%
%Note: This program needs to put following m_files in the same folder: 
%1. S (To define Skew Symmetric Matrix)
%2. C (To define Coriolis Matrix)
%3. D (To define Coriolis Matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%modify%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constant parameter

dh=0.01;                   % sampling time
t=0;                      % initial time
mult=6;                   % dimemsions of states 
multearth=6;
multtau=6;
kt=0;                     %index of array
x=[0;0;0;0;0;0];                  %SETTING INITIAL VALUE
xdt=[0;0;0;0;0;0];
eta=[0;0;0;0;0;0];
etadt=[0;0;0;0;0;0];
tau=[0;0;0;0;0;0];
taudt=[0;0;0;0;0;0];
    %%% Build M Matrix (Teetatbook Page 26.)
    m =35.6*10^6;  
    Zg = 0.0196;   
    Xg = 0;
    Yg=  0;
    rG=[Xg;Yg;Zg ];
    Zb = 0;   
    Xb = -0.611;
    Yb=  0;
    W=299;
    B=308;
           
    Ixx = 3.4*10^12;
    Izz = 60*10^12;
    Iyy = 60*10^12;
    Ixy=0;
    Ixz=0;
    Iyx=0;
    Iyz=0;
    Izx=0;
    Izy=0;
    I_0=[ Ixx  -Ixy  -Ixz ;
            -Iyx    Iyy  -Iyz ;  
            -Izx   -Izy   Izz ];
    
    M_11=m*eye(3);
    M_12=-m*S(rG);
    M_21=m*S(rG);
    M_22=I_0;
    M=[M_11 M_12; M_21 M_22];
    
%%%%%%%%%%%%%%%%%%%%%%%DO NOT CHANGE THIS PART%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ch=[.5 0. .5 0.];
cq=[2. 1. 1. 2.];
ckq=[.5 .29289322 1.7071068 .16666666];
ck=[.5 .29289322 1.7071068 .5];
qq=zeros(mult,1);           % relative to mult (dimemsion must be the same with state)
qq2=zeros(multearth,1);           % relative to mult (dimemsion must be the same with state)
qq3=zeros(multtau,1);           % relative to mult (dimemsion must be the same with state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%COMPUTING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while t<90
%% Define Time Varient Parameter
    %%% Build C Matrieta (Textbook Page 28)
    Cv=C(M,x);
    
    %%%Define D Matrieta(Textbook Page 43)
    Dv=D(x); % temporarily suppose
    
    %%%Define G Matrieta(Textbook Page 43)
    G_eta=G(eta); % temporarily suppose
    
    %%%Build J matrix
    C_etaphi=[ 1 , 0 , 0 ; 0 , cos(eta(4)) , sin(eta(4)) ; 0 , -sin(eta(4)) , cos(eta(4))];
    C_ytheta=[  cos(eta(5)) , 0 , -sin(eta(5)) ; 0 , 1 , 0 ; sin(eta(5)) , 0 , cos(eta(5))];
    C_zpsi=[  cos(eta(6)) , sin(eta(6)) , 0 ; -sin(eta(6)) , cos(eta(6)) , 0 ; 0 , 0 , 1];
    J_1=C_etaphi' * C_ytheta' * C_zpsi';

    J_2=[     1 , sin(eta(4))*tan(eta(5)) , cos(eta(4))*tan(eta(5)) ;
              0 ,          cos(eta(4))        ,           -sin(eta(4))      ;
              0 , sin(eta(4))/cos(eta(5)) , cos(eta(4))/cos(eta(5)) ];
    J=[J_1,zeros(3,3);zeros(3,3),J_2];
    
    %%%Build Jdot matrix     
    j11=-sin(eta(6))*etadt(6)*cos(eta(5))-sin(eta(5))*etadt(5)*cos(eta(6));
    j12=-cos(eta(6))*etadt(6)*cos(eta(4))+sin(eta(6))*sin(eta(4))*etadt(4)+(-sin(eta(6))*etadt(6)*sin(eta(5))+cos(eta(6))*sin(eta(5))*etadt(5))*sin(eta(4))+cos(eta(6))*sin(eta(5))*cos(eta(4))*etadt(4);
    j13= cos(eta(6))*etadt(6)*sin(eta(4))+sin(eta(6))*cos(eta(4))*etadt(4)+(-sin(eta(6))*etadt(6)*cos(eta(4))-cos(eta(6))*sin(eta(4))*etadt(4))*sin(eta(5))+cos(eta(6))*cos(eta(4))*cos(eta(5))*etadt(5);
    j21= cos(eta(6))*etadt(6)*cos(eta(5))-sin(eta(5))*etadt(5)*sin(eta(6));
    j22=-sin(eta(6))*etadt(6)*cos(eta(4))-sin(eta(4))*etadt(4)*cos(eta(6))+(cos(eta(4))*etadt(4)*sin(eta(5))+cos(eta(5))*etadt(5)*sin(eta(4)))*sin(eta(6))+sin(eta(4))*sin(eta(5))*cos(eta(6))*etadt(6);
    j23= sin(eta(6))*etadt(6)*sin(eta(4))-cos(eta(4))*etadt(4)*cos(eta(6))+(cos(eta(5))*etadt(5)*sin(eta(6))+cos(eta(6))*etadt(6)*sin(eta(5)))*cos(eta(4))-sin(eta(5))*sin(eta(6))*sin(eta(4))*etadt(4);
    j31=-cos(eta(5))*etadt(5);
    j32=-sin(eta(5))*etadt(5)*sin(eta(4))+cos(eta(4))*etadt(4)*cos(eta(5));
    j33=-sin(eta(5))*etadt(5)*cos(eta(4))-sin(eta(4))*etadt(4)*cos(eta(5));
    j45= cos(eta(4))*etadt(4)*tan(eta(5))+(sec(eta(5))^2)*etadt(5)*sin(eta(4));
    j46=-sin(eta(4))*etadt(4)*tan(eta(5))+(sec(eta(5))^2)*etadt(5)*cos(eta(4));
    j55=-sin(eta(4))*etadt(4);
    j56=-cos(eta(4))*etadt(4);
    j65= cos(eta(4))*etadt(4)*(cos(eta(5))^-1)+sin(eta(5))*etadt(5)*(cos(eta(5))^-2)*sin(eta(4));
    j66=-sin(eta(4))*etadt(4)*(cos(eta(5))^-1)+sin(eta(5))*etadt(5)*(cos(eta(5))^-2)*cos(eta(4));
  
    Jdot=[j11 j12 j13 0 0 0;j21 j22 j23 0 0 0; j31 j32 j33 0 0 0;0 0 0 0 j45 j46;0 0 0 0 j55 j56;0 0 0 0 j65 j66];
     
    for i=1:4             %System model start    
        %% Body-Frame State Space Equation (Mv' + C(v)v + D(v)v + g(£b) = £n)

            xdt=M\(tau-Cv*x-Dv*x-G_eta);

        %% Rotation Principle (Body -> Earth) (Textbook Page 48)

            etadt=J*x;
            etaddt= J*xdt+Jdot*x;
    
        %% Earth-Fietaed Frame Caculation
        
        
         %Control Input Actuator Time Delay 
%             A_thr=diag([-1/3,-1/3,-1/3,-1/3,-1/3,-1/3]);
            tau=[0;0;0;0;0;0];
%             taudt=A_thr*(tau-taud);
        %% Apply Current Force
             v_c=[0;0;0;0;0;0]; % current velocity
             v_r=x-v_c;
             Cvr=C(M,v_r);
             Dvr=D(v_r);
             tauc=M*xdt+Cvr*v_r+Dvr*v_r+G_eta;
             tau=tau+tauc;
               
        for mm=1:mult
             ak=dh*xdt(mm);
             rm=(ak-cq(i)*qq(mm))*ckq(i);
             x(mm)=x(mm)+rm;
             qq(mm)=qq(mm)+3.*rm-ck(i)*ak;   
        end
        
        for mm=1:multearth 
            ak2=dh*etadt(mm);
            rm2=(ak2-cq(i)*qq2(mm))*ckq(i);
            eta(mm)=eta(mm)+rm2;
            qq2(mm)=qq2(mm)+3.*rm2-ck(i)*ak2; 
        end
        for mm=1:multtau 
            ak3=dh*taudt(mm);
            rm3=(ak3-cq(i)*qq3(mm))*ckq(i);
            tau(mm)=tau(mm)+rm3;
            qq3(mm)=qq3(mm)+3.*rm3-ck(i)*ak3; 
        end
     t=t+dh*ch(i);
    end

kt=kt+1;
time(kt)=t;
ob(kt,:)=eta;
obbody(kt,:)=x;
obtau(kt,:)=tau;
end

%% Plot
figure(1)
plot3(ob(:,1),-ob(:,2),-ob(:,3))
figure(2)
plot(time,obbody(:,6))
figure(3)
plot(time,obtau(:,1))