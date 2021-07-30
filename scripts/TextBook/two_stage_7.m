%==================================================
%=     T W O - S T A G E  B O N D   G R A P H     =
%==================================================

%{
    .This team desires to create a transient model
        via the BOND GRAPH method
    .Eq(1)  qdot2   = Qin
    .Eq(2)  qdot5   = Qin - Qrad
    .Eq(3)  p11     = (Qrad + Qconv)/b11
    .Eq(4)  qdot2   = Qin + Qconv + Qrad - p11*b11

    .Current setup is to cook an 8oz steak no 8in cookplate

Author=         Joby Anthony III
Project=        Capstone 18 (19-20)
Date Created=   200204
%}

%{ 
REVISION & NOTES
    200203 - Dr. Medina suggested that we split the model to a two-stage model
    200221 - Josh and Joby talked through how to separate the processes
    200224 - Better separated time
    200225 - Re-meshed methods together
    200226 - Basically re did the whole thing
    200226(1) - Redid stage 2
    200302 - Learned that I was dumb and over-added the passage of time
    200324 - Edited the graphical output of results
%}


%% Clear Cache
clear,clc



%% Iterative Values
%Thermal mass
%Table Salt
rho_mass    = 135.469;                  %density of water [lb/ft3]
v_mass      = 0.0926;                   %volume of pyrex container [ft3]
m_mass      = rho_mass*v_mass;          %total mass of water [lb]
cp_mass     = 0.21;                     %specific heat of water [BTU/lb-F]

%Plate Radiation
%https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
epsilon     = 0.09;                     %emissivity of aluminum sheet

%Plate Conduction
%https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
k_plate     = 137;                      %thermal conductivity of aluminum [BTU/hr-ft-F]

%Peltier Circuit
%textbook value
k_pyrex     = 0.581;                    %thermal conductivity of pyrex [BTU/hr-ft-F]
%time of 74 min to approach operating limit w/o insulation
%https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
k_insul     = 0.0231;                   %thermal conductivity of fiberglass [BTU/hr-ft-F]
%no change in time w/ insulation
k_bot       = k_pyrex + k_insul;

%Conduction to Steak
moz         = 8;                        %weight of steak [oz]
%https://www.engineeringtoolbox.com/specific-heat-capacity-food-d_295.html
cp_steak    = 0.66;                     %specific heat of beef loin [BTU/lb-F]
%https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
k_steak     = 0.248;                    %thermal conductivity of beef loin [BTU/hr-ft-F]



%% Constants and Parameters
%Initial Temperatures
TsurfI      = 60;                       %initialize at T_not [F]
TsprevI     = TsurfI;
TmassI      = TsurfI;
Tnot        = 60;
TinfI       = 72;                       %bulk fluid temperature [F]
TsteakI     = 32;                       %initial temperature of steak [F]
TtprevI     = TsteakI;

%Input Flux from Sun
%https://www.sciencedirect.com/topics/earth-and-planetary-sciences/solar-flux
%yes, joby. that "3600" is supposed to be there.
%STOP DELETING THIS LINE!!!
Qsun        = 0.1209*3600;              %heat flux of sun for work day [Btu/hr-ft2]
Afl         = 1.2917;                   %area of Fresnel Lens [ft2]
t           = 0.92;                     %transmittance [%]

%Thermal Mass Constant
b_mass      = m_mass*cp_mass;           %capacitance constant of mass [BTU/F]

%Steak Constant
mlb         = moz/16;                   %weight of steak [lb]
b_steak     = mlb*cp_steak;             %capacitance constant of steak [BTU/F]
A_steak     = 0.5;                      %area of steak [ft2]
x_steak     = 1/12;                     %thickness of steak [ft]
rk_steak    = k_steak*A_steak/x_steak;  %resistor constant of steak [BTU/hr-F]

%Qcond Constant
A_plate     = 0.5625;                   %initial area of cook plate [ft2]
x_plate     = 0.25/12;                  %thickness of cook plate [ft]
rk_plate    = k_plate*A_plate/x_plate;  %resistor constant of plate [BTU/hr-F]

%Qconv Modulushc
UwindI      = 5.2;%:0.1:7.8;            %cross-wind velocity [mph]
MphToMs     = 0.44704;                  %convert mph to m/s
UwindK      = UwindI.*MphToMs;          %cross-wind velocity [m/s]
WToBTUs     = 1/1055.055852;            %W to BTU/s
hc          = 1.16.*(10.45 - UwindK...
                + 10.*sqrt(UwindK)).*0.000048919;    %film coefficient [BTU/hr-ft2-F]
m_conv      = hc*A_plate;               %transfer modulus [BTU/hr-F]

%Radiation Modulus
sb          = 1.714e-9;                 %stefan-boltzmann constant [BTU/hr-ft2-R4]
m_rad       = epsilon*sb*A_plate;       %radiation modulus [BTU/hr-R4]

%qdot5
dx          = 0.25/12;                  %thickness of bottom of Pyrex [ft]
T52         = 302;                      %operating limit of Peltier Circuit [F]
rk_peltier  = k_bot*A_plate/dx;



%% Initialize Calculations
%Qin
Qin         = t*Qsun*Afl;               %heat flux input [BTU/hr]

time        = 2;                        %time range of analysis [hr]
n           = 100;                     %number of desired iterations
k           = time/n;                   %time between iterations

Tpeltier = 60;
%T_peltier   = zeros(1,n);
%TMASSI      = zeros(1,n);
%t_mass      = zeros(1,n);
%TSURFI      = zeros(1,n);

%TSTEAKI     = zeros(1,n);
%TMASSI_prime= zeros(1,n);
%t_steak     = zeros(1,n);
%TSURFI_prime= zeros(1,n);



%% Thermal Mass Temperature
%Solve for plate surface temperature
fprintf('STAGE 1\n')
fprintf('Iteration  Mass     Surface    Peltier     Time\n')
fprintf('-----------------------------------------------\n')
a = 0;
for i = 1:1:n
    if Tpeltier < T52
    
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Recall and convert temperatures
        TmprevI     = TmassI;
        TsprevI     = TsurfI;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve temperature and time of current iteration
        t_mass(i)       = (i-1)*k*60;         %array of total time [min]
        TmassI      = Qin*k/b_mass + TmprevI; %[F]
        TMASSI(i)       = TmprevI;    


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve for the mass temperature with net loss heat flux
        % Qnet = Qin - (Qrad + Qconv)
        Qtm         = b_mass*(TmassI - TmprevI);
        Tpeltier= Qtm/rk_peltier + TmprevI;
        T_peltier(i) = Tpeltier;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve the new plate temperature
        rhs = Qtm;                          %heat flux into system [BTU]
        syms srk_plate sm_rad sm_conv sTmassI sTsprevI sTinfI sk
        % Qtm = Qcond - (Qrad + Qconv)
        eqn = (srk_plate*(sTmassI - sTsprevI) - ...
            (sm_rad*((sTsprevI + 459.67)^4 - (sTinfI + 459.67)^4) + ...
            sm_conv*(sTsprevI - sTinfI)))*sk == rhs;
        s = solve(eqn, sTsprevI, 'MaxDegree', 4);
        old = [srk_plate sm_rad sm_conv sTmassI sTsprevI sTinfI sk];
        new = [rk_plate m_rad m_conv TmassI TsprevI TinfI k];
        arraySolve = abs(double(subs(s,old,new)));
        TsurfI = min(arraySolve);
        TSURFI(i)   = TsprevI;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Print number of iterations performed
        fprintf ('%9.0f     %4.1f   %7.1f   %7.1f   %4.2f\n',...
            i-1,TmprevI,TsprevI,Tpeltier,t_mass(i));
        a = a + 1;
    end
end

%TMASSI = size(1,a);
%TSURFI = size(1,a);
%T_peltier = size(1,a);

fprintf('\n')
f = figure(1);
plot(t_mass,TMASSI,'r-','DisplayName','Thermal Mass')
hold on
plot(t_mass,TSURFI,'g-','DisplayName','Cook Plate')
plot(t_mass,T_peltier,'b-','DisplayName','Peltier Circuit')
grid on
title 'Temperature of Thermal Mass (Water)'
xlabel 'Time [min]'
ylabel 'Temperature [°F]'
legend('show','Location','northwest')
hold off



%% Surface and Steak Temperature
Q_prime     = b_mass*(TmprevI - Tnot);

%Solve for steak temperature
TmassI_prime= TmprevI;
TsurfI_prime= TsprevI;

fprintf('STAGE 2\n')
fprintf('Iteration  Mass     Surface    Steak   Time\n')
fprintf('-------------------------------------------\n')
b = 0;
for j = 1:1:n
    if TmassI_prime > Tnot
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Initialize steak and recall by iteration
        TtprevI     = TsteakI;
        TmprevI_prime=TmassI_prime;
        TsprevI_prime=TsurfI_prime;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve heat transfer of current iteration
        t_steak(j)      = (j-1)*k*60;



        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve for Qmass
        Qcond_plate = rk_plate*k*(TsprevI_prime - TmprevI_prime);
        %Qtm_prime   = Qcond_plate;
        TmassI_prime= -Q_prime*k/b_mass + TmprevI_prime;
        TMASSI_prime(j) = TmprevI_prime;
        %Qsteak      = rk_steak*tsteak*(TtprevI - TsprevI_prime);
        %Qnet        = Qcond_plate - ...
        %    (m_rad*((TsprevI_prime + 459.67)^4 - (TinfI + 459.67)^4) + ...
        %    m_conv*(TsprevI_prime - TinfI))*tsteak;
        TsteakI     = Q_prime*k/rk_steak + TtprevI;
        QPRIME(j)       = Q_prime;
        TSTEAKI(j)  = TtprevI;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Solve the new plate temperature
        rhs=Qcond_plate; %heat flux into system [B
        syms srk_plate sm_rad sm_conv sTmprevI_prime sTsprevI_prime sTtprevI sTinfI sk
        eqn = (srk_plate*(sTmprevI_prime - sTsprevI_prime) - ...
            (sm_rad*((sTsprevI_prime + 459.67)^4 - (sTinfI + 459.67)^4) + ...
            sm_conv*(sTsprevI_prime - sTinfI) + (sTtprevI - sTinfI)))*sk == rhs;
        s = solve(eqn, sTsprevI_prime, 'MaxDegree', 4);
        old = [srk_plate sm_rad sm_conv sTmprevI_prime sTsprevI_prime sTtprevI sTinfI sk];
        new = [rk_plate m_rad m_conv TmprevI_prime TsprevI_prime TtprevI TinfI k];
        arraySolve = abs(double(subs(s,old,new)));
        TsurfI_prime = min(arraySolve);
        TSURFI_prime(j) = TsprevI_prime;


        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        %Print number of iterations performed
        fprintf ('%9.0f     %4.1f   %7.1f   %6.1f   %4.2f\n',...
            j-1,TmprevI_prime,TsprevI_prime,TtprevI,t_steak(j));
        b = b + 1;
    end
end

%TSTEAKI = size(1,a);
%TSURFI_prime = size(1,a);
%TMASSI_prime = size(1,a);

g = figure(2);
plot(t_steak,TSTEAKI,'g-','DisplayName','Food Item')
hold on
plot(t_steak,TSURFI_prime,'k-','DisplayName','Cook Plate')
plot(t_steak,TMASSI_prime,'b-','DisplayName','Thermal Mass')
grid on
title 'Temperature of Food Item (Steak)'
xlabel 'Time [min]' 
ylabel 'Temperature [°F]'
legend('show')
hold off



%% End of Code
%end of code


