%% Defining Transfer Functions for Conventional and Third Damper Systems

% Define symbolic variables
syms m_us1 m_s m_us2 k_t1 k_t2 k1 k b b1 b2 k2 s omega;

%==============Conventional Suspension System======================%
A = [   0,   0,   0,     0,  0,  -1/m_us1,           0,         0; 
        0,   0,   0,     0,  0,   1/m_us1,      -1/m_s,         0; 
        0,   0,   0,     0,  0,         0,      -1/m_s,   1/m_us2; 
        0,   0,   0,     0,  0,         0,           0,   1/m_us2; 
        0,   0,   0,     0,  0,   1/m_us1,           0,  -1/m_us2; 
     k_t1, -k1,   0,     0, -k, -b1/m_us1,      b1/m_s,         0; 
        0,  k1,  k2,     0,  0,  b1/m_us1,(-b1-b2)/m_s,  b2/m_us2; 
        0,   0, -k2, -k_t2,  k,         0,      b2/m_s, -b2/m_us2];
B = [1, 0, 0, 0; 
     0, 0, 0, 0; 
     0, 0, 0, 0; 
     0, 0, 0, 0; 
     0, 0, 0, 0; 
     0, -1, 0, 0; 
     0, 0, 0, -1; 
     0, 0, -1, 0];

A_s = s*eye(8)-A;

% Transfer Function for Q2/Vin
Num1 = [B(:,1),A_s(:,2:8)];
Q2_Vin = det(Num1)/det(A_s);

% Transfer Function for Q8/Vin
Num2 = [A_s(:,1),B(:,1),A_s(:,3:8)];
Q8_Vin = det(Num2)/det(A_s);

% Transfer Function for Q15/Vin
Num3 = [A_s(:,1:2),B(:,1),A_s(:,4:8)];
Q15_Vin = det(Num3)/det(A_s);

% Transfer Function for Q20/Vin
Num4 = [A_s(:,1:3),B(:,1),A_s(:,5:8)];
Q20_Vin = det(Num4)/det(A_s);

% Transfer Function for Q22/Vin
Num5 = [A_s(:,1:4),B(:,1),A_s(:,6:8)];
Q22_Vin = det(Num5)/det(A_s);

% Transfer Function for P4/Vin
Num6 = [A_s(:,1:5),B(:,1),A_s(:,7:8)];
P4_Vin = det(Num6)/det(A_s);

% Transfer Function for P12/Vin
Num7 = [A_s(:,1:6),B(:,1),A_s(:,8)];
P12_Vin = det(Num7)/det(A_s);

% Transfer Function for P19/Vin
Num8 = [A_s(:,1:7),B(:,1)];
P19_Vin = det(Num8)/det(A_s);

% Transfer Function for F4/Vin
F4plusm_us1g_Vin = k_t1*Q2_Vin - k1*Q8_Vin - k*Q22_Vin - b1/m_us1*P4_Vin + b1/m_s*P12_Vin;

% Transfer Function for F12/Vin
F12plusm_sg_Vin = k1*Q8_Vin + k2*Q15_Vin + b1/m_us1*P4_Vin + (-b1-b2)/m_s*P12_Vin + b2/m_us2*P19_Vin;

% Transfer Function for F19_Vin
F19plusm_us2g_Vin = -k2*Q15_Vin - k_t2*Q20_Vin + k*Q22_Vin + b2/m_s*P12_Vin - b2/m_us2*P19_Vin;

% Substitute j*omega for s and make dimensionless by dividing by input tire
% damping coefficient
F4plusm_us1g_Vinb1_O = subs(F4plusm_us1g_Vin, s, omega*1i)/b1;
F12plusm_sg_Vinb1_O = subs(F12plusm_sg_Vin, s, omega*1i)/b1;
F19plusm_us2g_Vinb1_O = subs(F19plusm_us2g_Vin, s, omega*1i)/b1;

%==============Transverse Third Damper System===================%
A_t = [   0,   0,   0,     0,  0,      -1/m_us1,            0,              0; 
          0,   0,   0,     0,  0,       1/m_us1,       -1/m_s,              0; 
          0,   0,   0,     0,  0,             0,       -1/m_s,        1/m_us2; 
          0,   0,   0,     0,  0,             0,            0,        1/m_us2; 
          0,   0,   0,     0,  0,       1/m_us1,            0,       -1/m_us2; 
       k_t1, -k1,   0,     0, -k, (-b1-b)/m_us1,       b1/m_s,       -b/m_us2; 
          0,  k1,  k2,     0,  0,      b1/m_us1, (-b1-b2)/m_s,       b2/m_us2; 
          0,   0, -k2, -k_t2,  k,      -b/m_us1,       b2/m_s, (-b-b2)/m_us2];
B_t = [1, 0, 0, 0; 
       0, 0, 0, 0; 
       0, 0, 0, 0; 
       0, 0, 0, 0; 
       0, 0, 0, 0; 
       0, -1, 0, 0; 
       0, 0, 0, -1; 
       0, 0, -1, 0];

A_s_t = s*eye(8)-A_t;

% Transfer Function for Q2/Vin
Num1_t = [B_t(:,1),A_s_t(:,2:8)];
Q2_Vin_t = det(Num1_t)/det(A_s_t);

% Transfer Function for Q8/Vin
Num2_t = [A_s_t(:,1),B_t(:,1),A_s_t(:,3:8)];
Q8_Vin_t = det(Num2_t)/det(A_s_t);

% Transfer Function for Q15/Vin
Num3_t = [A_s_t(:,1:2),B_t(:,1),A_s_t(:,4:8)];
Q15_Vin_t = det(Num3_t)/det(A_s_t);

% Transfer Function for Q20/Vin
Num4_t = [A_s_t(:,1:3),B_t(:,1),A_s_t(:,5:8)];
Q20_Vin_t = det(Num4_t)/det(A_s_t);

% Transfer Function for Q22/Vin
Num5_t = [A_s_t(:,1:4),B_t(:,1),A_s_t(:,6:8)];
Q22_Vin_t = det(Num5_t)/det(A_s_t);

% Transfer Function for P4/Vin
Num6_t = [A_s_t(:,1:5),B_t(:,1),A_s_t(:,7:8)];
P4_Vin_t = det(Num6_t)/det(A_s_t);

% Transfer Function for P12/Vin
Num7_t = [A_s_t(:,1:6),B_t(:,1),A_s_t(:,8)];
P12_Vin_t = det(Num7_t)/det(A_s_t);

% Transfer Function for P19/Vin
Num8_t = [A_s_t(:,1:7),B_t(:,1)];
P19_Vin_t = det(Num8_t)/det(A_s_t);

% Transfer Function for F4/Vin
F4plusm_us1g_Vin_t = k_t1*Q2_Vin_t - k1*Q8_Vin_t - k*Q22_Vin_t - (b1+b)/m_us1*P4_Vin_t + b1/m_s*P12_Vin_t - b/m_us2*P19_Vin_t;

% Transfer Function for F12/Vin
F12plusm_sg_Vin_t = k1*Q8_Vin_t + k2*Q15_Vin_t + b1/m_us1*P4_Vin_t + (-b1-b2)/m_s*P12_Vin_t + b2/m_us2*P19_Vin_t;

% Transfer Function for F19_Vin
F19plusm_us2g_Vin_t = -k2*Q15_Vin_t - k_t2*Q20_Vin_t + k*Q22_Vin_t - b/m_us1*P4_Vin_t + b2/m_s*P12_Vin_t - (b2+b)/m_us2*P19_Vin_t;

% Substitute j*omega for s and make dimensionless by dividing by input tire
% damping coefficient
F4plusm_us1g_Vinb1_t_O = subs(F4plusm_us1g_Vin_t, s, omega*1i)/b1;
F12plusm_sg_Vinb1_t_O = subs(F12plusm_sg_Vin_t, s, omega*1i)/b1;
F19plusm_us2g_Vinb1_t_O = subs(F19plusm_us2g_Vin_t, s, omega*1i)/b1;

%% Define electric car parameters

M_s = 680; %kg
M_us1 = 25; %kg
M_us2 = M_us1; %kg
K1 = 32000; %N/m
K2 = K1; %N/m
K_t1 = 360000; %N/m
K_t2 = K_t1; %N/m
K = 360000; %N/m
DR = 0.6;
B1 = DR*2*sqrt((M_s/4)*K1);
B2 = B1;
C = 0.4;
B = C*B1;

%% Plot unsprung mass effect on AR of sprung mass

% Set unsprung mass to low
M_us1 = 25; %kg
M_us2 = M_us1;

% Calculate AR of sprung mass for conventional
F12plusm_sg_Vinb1 = subs(F12plusm_sg_Vinb1_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B1,B2]);
[F12plusm_sg_Vinb1_num,F12plusm_sg_Vinb1_den] = numden(F12plusm_sg_Vinb1);
amp_ratio_F12 = sqrt((real(F12plusm_sg_Vinb1_num))^2+(imag(F12plusm_sg_Vinb1_num))^2)/sqrt((real(F12plusm_sg_Vinb1_den))^2+(imag(F12plusm_sg_Vinb1_den))^2);

% Calculate AR of sprung mass for transverse 3rd damper
F12plusm_sg_Vinb1_t = subs(F12plusm_sg_Vinb1_t_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B,B1,B2]);
[F12plusm_sg_Vinb1_t_num,F12plusm_sg_Vinb1_t_den] = numden(F12plusm_sg_Vinb1_t);
amp_ratio_F12_t = sqrt((real(F12plusm_sg_Vinb1_t_num))^2+(imag(F12plusm_sg_Vinb1_t_num))^2)/sqrt((real(F12plusm_sg_Vinb1_t_den))^2+(imag(F12plusm_sg_Vinb1_t_den))^2);

% Calculate phase angle for conventional
phi_F12 = atan2d((imag(F12plusm_sg_Vinb1_num)),real(F12plusm_sg_Vinb1_den))-atan2d(imag(F12plusm_sg_Vinb1_den),real(F12plusm_sg_Vinb1_den));
phi_F12 = double(subs(phi_F12,omega,0:1000));
phi_F12a = phi_F12.*(phi_F12 >=0)+(phi_F12+360).*(phi_F12 <0);
for i = 11:length(phi_F12a)
    if phi_F12a(i) > 1
        phi_F12a(i) = phi_F12a(i) - 360;
    else
        phi_F12a(i) = phi_F12a(i);
    end
end

% Calculate phase angle for transverse 3rd damper
phi_F12_t = atan2d((imag(F12plusm_sg_Vinb1_t_num)),real(F12plusm_sg_Vinb1_t_num))-atan2d(imag(F12plusm_sg_Vinb1_t_den),real(F12plusm_sg_Vinb1_t_den));
phi_F12_t = double(subs(phi_F12_t,omega,0:1000));
phi_F12_t_a = phi_F12_t.*(phi_F12_t >=0)+(phi_F12_t+360).*(phi_F12_t <0);
for i = 12:length(phi_F12_t_a)
    if phi_F12_t_a(i) > 50
        phi_F12_t_a(i) = phi_F12_t_a(i) - 360;
    else
        phi_F12_t_a(i) = phi_F12_t_a(i);
    end
end

% Create the plot
subplot(2,1,1)
fplot(amp_ratio_F12,[1,1000],'-r');
title('Unsprung Mass Effects on AR of Sprung Mass','FontSize',18);
xlabel('Frequency [$Hz$]','Interpreter','LaTeX','FontSize',18);
ylabel('Amplitude Ratio [$\frac{F_{ms}+m_{s}g}{V_{in}b_1}$]','Interpreter','LaTeX','FontSize',18);
ylim([0,4])
hold on
fplot(amp_ratio_F12_t,[1,1000],'-b');

% Set unsprung mass to MEDIUM
M_us1 = 50; %kg
M_us2 = M_us1;

% Recalculate AR for conventional
F12plusm_sg_Vinb1 = subs(F12plusm_sg_Vinb1_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B1,B2]);
[F12plusm_sg_Vinb1_num,F12plusm_sg_Vinb1_den] = numden(F12plusm_sg_Vinb1);
amp_ratio_F12 = sqrt((real(F12plusm_sg_Vinb1_num))^2+(imag(F12plusm_sg_Vinb1_num))^2)/sqrt((real(F12plusm_sg_Vinb1_den))^2+(imag(F12plusm_sg_Vinb1_den))^2);

% Recalculate AR for transverse 3rd damper
F12plusm_sg_Vinb1_t = subs(F12plusm_sg_Vinb1_t_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B,B1,B2]);
[F12plusm_sg_Vinb1_t_num,F12plusm_sg_Vinb1_t_den] = numden(F12plusm_sg_Vinb1_t);
amp_ratio_F12_t = sqrt((real(F12plusm_sg_Vinb1_t_num))^2+(imag(F12plusm_sg_Vinb1_t_num))^2)/sqrt((real(F12plusm_sg_Vinb1_t_den))^2+(imag(F12plusm_sg_Vinb1_t_den))^2);

% Recalculate phase angle for conventional
phi_F12 = atan2d((imag(F12plusm_sg_Vinb1_num)),real(F12plusm_sg_Vinb1_den))-atan2d(imag(F12plusm_sg_Vinb1_den),real(F12plusm_sg_Vinb1_den));
phi_F12 = double(subs(phi_F12,omega,0:1000));
phi_F12b = phi_F12.*(phi_F12 >=0)+(phi_F12+360).*(phi_F12 <0);
for i = 11:length(phi_F12b)
    if phi_F12b(i) > 1
        phi_F12b(i) = phi_F12b(i) - 360;
    else
        phi_F12b(i) = phi_F12b(i);
    end
end

% Calculate phase angle for transverse 3rd damper
phi_F12_t = atan2d((imag(F12plusm_sg_Vinb1_t_num)),real(F12plusm_sg_Vinb1_t_num))-atan2d(imag(F12plusm_sg_Vinb1_t_den),real(F12plusm_sg_Vinb1_t_den));
phi_F12_t = double(subs(phi_F12_t,omega,0:1000));
phi_F12_t_b = phi_F12_t.*(phi_F12_t >=0)+(phi_F12_t+360).*(phi_F12_t <0);
for i = 12:length(phi_F12_t_b)
    if phi_F12_t_b(i) > 50
        phi_F12_t_b(i) = phi_F12_t_b(i) - 360;
    else
        phi_F12_t_b(i) = phi_F12_t_b(i);
    end
end

% Create plot
fplot(amp_ratio_F12,[1,1000],'--r');
fplot(amp_ratio_F12_t,[1,1000],'--b');

% Set unsprung mass to HIGH
M_us1 = 75; %kg
M_us2 = M_us1;

% Recalculate AR for conventional
F12plusm_sg_Vinb1 = subs(F12plusm_sg_Vinb1_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B1,B2]);
[F12plusm_sg_Vinb1_num,F12plusm_sg_Vinb1_den] = numden(F12plusm_sg_Vinb1);
amp_ratio_F12 = sqrt((real(F12plusm_sg_Vinb1_num))^2+(imag(F12plusm_sg_Vinb1_num))^2)/sqrt((real(F12plusm_sg_Vinb1_den))^2+(imag(F12plusm_sg_Vinb1_den))^2);

% Recalculate AR for transverse 3rd damper
F12plusm_sg_Vinb1_t = subs(F12plusm_sg_Vinb1_t_O, [m_s,m_us1,m_us2,k1,k2,k_t1,k_t2,k,b,b1,b2], [M_s,M_us1,M_us2,K1,K2,K_t1,K_t2,K,B,B1,B2]);
[F12plusm_sg_Vinb1_t_num,F12plusm_sg_Vinb1_t_den] = numden(F12plusm_sg_Vinb1_t);
amp_ratio_F12_t = sqrt((real(F12plusm_sg_Vinb1_t_num))^2+(imag(F12plusm_sg_Vinb1_t_num))^2)/sqrt((real(F12plusm_sg_Vinb1_t_den))^2+(imag(F12plusm_sg_Vinb1_t_den))^2);

% Recalculate phase angle for conventional
phi_F12 = atan2d((imag(F12plusm_sg_Vinb1_num)),real(F12plusm_sg_Vinb1_den))-atan2d(imag(F12plusm_sg_Vinb1_den),real(F12plusm_sg_Vinb1_den));
phi_F12 = double(subs(phi_F12,omega,0:1000));
phi_F12c = phi_F12.*(phi_F12 >=0)+(phi_F12+360).*(phi_F12 <0);
for i = 11:length(phi_F12c)
    if phi_F12c(i) > 1
        phi_F12c(i) = phi_F12c(i) - 360;
    else
        phi_F12c(i) = phi_F12c(i);
    end
end

% Recalculate phase angle for transverse 3rd damper
phi_F12_t = atan2d((imag(F12plusm_sg_Vinb1_t_num)),real(F12plusm_sg_Vinb1_t_num))-atan2d(imag(F12plusm_sg_Vinb1_t_den),real(F12plusm_sg_Vinb1_t_den));
phi_F12_t = double(subs(phi_F12_t,omega,0:1000));
phi_F12_t_c = phi_F12_t.*(phi_F12_t >=0)+(phi_F12_t+360).*(phi_F12_t <0);
for i = 12:length(phi_F12_t_c)
    if phi_F12_t_c(i) > 50
        phi_F12_t_c(i) = phi_F12_t_c(i) - 360;
    else
        phi_F12_t_c(i) = phi_F12_t_c(i);
    end
end

% Create plot
fplot(amp_ratio_F12,[0,1000],':r');
fplot(amp_ratio_F12_t,[0,1000],':b');
hold off
legend({'Conventional: Low Unsprung Mass (25kg)','3rd Damper: Low Unsprung Mass (25 kg)','Conventional: Medium Unsprung Mass (50 kg)','3rd Damper: Medium Unsprung Mass (50 kg)','Conventional: High Unsprung Mass (75 kg)','3rd Damper: High Unsprung Mass (75 kg)'},'Location','northwest','FontSize',10);
set(gca,'XScale','log');

% Plot phase angles
subplot (2,1,2);
plot(0:1000,phi_F12a,'-r');
title('Unsprung Mass Effects on Phase Angle of Sprung Mass','FontSize',15);
xlabel('Frequency [$Hz$]','Interpreter','LaTeX','FontSize',17);
ylabel('Phase Angle [$^{\circ}$]','Interpreter','LaTeX','FontSize',17);
hold on
plot(0:1000,phi_F12_t_a,'-b');
plot(0:1000,phi_F12b,'--r');
plot(0:1000,phi_F12_t_b,'--b');
plot(0:1000,phi_F12c,':r');
plot(0:1000,phi_F12_t_c,':b');
hold off
set(gca,'XScale','log');
