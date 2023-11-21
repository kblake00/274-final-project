clear
name = 'walker';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t th1 th2 dth1 dth2 ddth1 ddth2 th3 th4 dth3 dth4 ddth3 ddth4 x y dx dy ddx ddy real
syms m1 m2 m3 m4 I1 I2 I3 I4 l_O_m1 l_B_m2 l_A_m3 l_C_m4 g ramp real
syms l_OA l_OB l_AC l_DE real 
syms tau1 tau2 tau3 tau4 Fx Fy real
syms Ir N real

% Group them
q   = [x; y; th1  ; th2 ; th3; th4];      % generalized coordinates
dq  = [dx; dy; dth1 ; dth2; dth3; dth4];    % first time derivatives
ddq = [ddx; ddy; ddth1;ddth2; ddth3; ddth4];  % second time derivatives
u   = [tau1 ; tau2; tau3; tau4];     % controls
F   = [Fx ; Fy];

p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g ramp]';        % parameters

% Generate Vectors and Derivatives
ihat = [0; -1; 0];
jhat = [1; 0; 0];

xhat = [1; 0; 0];
yhat = [0; 1; 0];

iramp = jhat*cos(ramp) -ihat*sin(ramp);
jramp = -jhat*sin(ramp) - ihat*cos(ramp);

khat = cross(ihat,jhat);
e1hat_l =  cos(th1)*ihat + sin(th1)*jhat;
e2hat_l =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;

e1hat_r =  cos(th3)*ihat + sin(th3)*jhat;
e2hat_r =  cos(th3+th4)*ihat + sin(th3+th4)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

r0 = x*jhat - y*ihat;

rA_l = r0 + l_OA * e1hat_l;
rB_l = r0 + l_OB * e1hat_l;
rC_l = rA_l  + l_AC * e2hat_l;
rD_l = rB_l  + l_AC * e2hat_l;
rE_l = rD_l  + l_DE * e1hat_l;

r_m1_l = r0 + l_O_m1 * e1hat_l;
r_m2_l = rB_l + l_B_m2 * e2hat_l;
r_m3_l = rA_l + l_A_m3 * e2hat_l;
r_m4_l = rC_l + l_C_m4 * e1hat_l;

dr0   = ddt(r0);
%ddr0 = simplify(ddt(dr0));

drA_l = ddt(rA_l);
drB_l = ddt(rB_l);
drC_l = ddt(rC_l);
drD_l = ddt(rD_l);
drE_l = ddt(rE_l);

ddrE_l = ddt(drE_l);

dr_m1_l = ddt(r_m1_l);
dr_m2_l = ddt(r_m2_l);
dr_m3_l = ddt(r_m3_l);
dr_m4_l = ddt(r_m4_l);

% now repeat everything for right leg
rA_r = r0 + l_OA * e1hat_r;
rB_r = r0 + l_OB * e1hat_r;
rC_r = rA_r  + l_AC * e2hat_r;
rD_r = rB_r  + l_AC * e2hat_r;
rE_r = rD_r  + l_DE * e1hat_r;

r_m1_r = r0 + l_O_m1 * e1hat_r;
r_m2_r = rB_r + l_B_m2 * e2hat_r;
r_m3_r = rA_r + l_A_m3 * e2hat_r;
r_m4_r = rC_r + l_C_m4 * e1hat_r;

drA_r = ddt(rA_r);
drB_r = ddt(rB_r);
drC_r = ddt(rC_r);
drD_r = ddt(rD_r);
drE_r = ddt(rE_r);

ddrE_r = ddt(drE_r);

dr_m1_r = ddt(r_m1_r);
dr_m2_r = ddt(r_m2_r);
dr_m3_r = ddt(r_m3_r);
dr_m4_r = ddt(r_m4_r);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces


omega1_l = dth1;
omega2_l = dth1 + dth2;
omega3_l = dth1 + dth2;
omega4_l = dth1;

T0 = (1/2)*m1 * dot(dr0,dr0); % ARBITRARILY BODY HAS SAME MASS AS M1

T1_l = (1/2)*m1 * dot(dr_m1_l,dr_m1_l) + (1/2) * I1 * omega1_l^2;
T2_l = (1/2)*m2 * dot(dr_m2_l,dr_m2_l) + (1/2) * I2 * omega2_l^2;
T3_l = (1/2)*m3 * dot(dr_m3_l,dr_m3_l) + (1/2) * I3 * omega3_l^2;
T4_l = (1/2)*m4 * dot(dr_m4_l,dr_m4_l) + (1/2) * I4 * omega4_l^2;
T1r_l = (1/2)*Ir*(N*dth1)^2;
T2r_l = (1/2)*Ir*(dth1 + N*dth2)^2;

Vg1_l = m1*g*dot(r_m1_l, jramp);
Vg2_l = m2*g*dot(r_m2_l, jramp);
Vg3_l = m3*g*dot(r_m3_l, jramp);
Vg4_l = m4*g*dot(r_m4_l, jramp);

% now repeat everything for right leg
omega1_r = dth3;
omega2_r = dth3 + dth4;
omega3_r = dth3 + dth4;
omega4_r = dth3;

T1_r = (1/2)*m1 * dot(dr_m1_r,dr_m1_r) + (1/2) * I1 * omega1_r^2;
T2_r = (1/2)*m2 * dot(dr_m2_r,dr_m2_r) + (1/2) * I2 * omega2_r^2;
T3_r = (1/2)*m3 * dot(dr_m3_r,dr_m3_r) + (1/2) * I3 * omega3_r^2;
T4_r = (1/2)*m4 * dot(dr_m4_r,dr_m4_r) + (1/2) * I4 * omega4_r^2;
T1r_r = (1/2)*Ir*(N*dth3)^2;
T2r_r = (1/2)*Ir*(dth3 + N*dth4)^2;

Vg0 = m1*g*dot(r0, jramp); % ARBITRARILY BODY HAS SAME MASS AS M1
Vg1_r = m1*g*dot(r_m1_r, jramp);
Vg2_r = m2*g*dot(r_m2_r, jramp);
Vg3_r = m3*g*dot(r_m3_r, jramp);
Vg4_r = m4*g*dot(r_m4_r, jramp);

T = simplify(T0 + T1_l + T2_l + T3_l + T4_l + T1r_l + T2r_l + T1_r + T2_r + T3_r + T4_r + T1r_r + T2r_r);
Vg = Vg0 + Vg1_l + Vg2_l + Vg3_l + Vg4_l + Vg1_r + Vg2_r + Vg3_r + Vg4_r;

% tau1 --> th1, tau2 --> th2 (left); tau3 --> th3, tau4 --> th4 (right); 
Q_tau1_l = M2Q(tau1*khat,omega1_l*khat);
Q_tau2_l = M2Q(tau2*khat,omega2_l*khat); 
Q_tau2R_l= M2Q(-tau2*khat,omega1_l*khat);

Q_tau1_r = M2Q(tau3*khat,omega1_r*khat);
Q_tau2_r = M2Q(tau4*khat,omega2_r*khat); 
Q_tau2R_r= M2Q(-tau4*khat,omega1_r*khat);


Q_tau = Q_tau1_l+Q_tau2_l + Q_tau2R_l + Q_tau1_r + Q_tau2_r + Q_tau2R_r;

Q = Q_tau;

% Assemble the array of cartesian coordinates of the key points
keypoints = [r0(1:2) rA_l(1:2) rB_l(1:2) rC_l(1:2) rD_l(1:2) rE_l(1:2) rA_r(1:2) rB_r(1:2) rC_r(1:2) rD_r(1:2) rE_r(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+Vg;
L = T-Vg;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;


% Rearrange Equations of Motion
A = simplify(jacobian(eom,ddq));
b = simplify(A*ddq - eom);

% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(Vg, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);

% Compute foot jacobian
J_l = jacobian(rE_l,q);
J_r = jacobian(rE_r,q);

% Compute ddt( J )
dJ_l= reshape( ddt(J_l(:)) , size(J_l) );
dJ_r= reshape( ddt(J_r(:)) , size(J_r) );

%Lambda_l = pinv(J_l*(pinv(A)*transpose(J_l)));
%Mu_l = Lambda_l*(J_l*(pinv(A)*Corr_Joint_Sp)) - Lambda_l*(dJ_l*dq);
%Rho_l = Lambda_l*(J_l*(pinv(A)*Grav_Joint_Sp));

%Lambda_r = pinv(J_r*(pinv(A)*transpose(J_r)));
%Mu_r = Lambda_r*(J_r*(pinv(A)*Corr_Joint_Sp)) - Lambda_r*(dJ_r*dq);
%Rho_r = Lambda_r*(J_r*(pinv(A)*Grav_Joint_Sp));

% Write Energy Function and Equations of Motion
z  = [q ; dq];

rE_l = rE_l(1:2);
drE_l= drE_l(1:2);
J_l  = J_l(1:2,:);
dJ_l = dJ_l(1:2,:);

rE_r = rE_r(1:2);
drE_r= drE_r(1:2);
J_r  = J_r(1:2,:);
dJ_r = dJ_r(1:2,:);

matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});

% Body info stuff
matlabFunction(r0,'file',['r0_' name],'vars',{z p});
matlabFunction(dr0,'file',['dr0_' name],'vars',{z p});
%matlabFunction(ddr0,'file',['ddr0_' name],'vars',{z p});



matlabFunction(rE_l,'file',['position_foot_left'],'vars',{z p});
matlabFunction(drE_l,'file',['velocity_foot_left'],'vars',{z p});
%matlabFunction(ddrE_l,'file',['acceleration_foot_left'],'vars',{z p});
matlabFunction(J_l ,'file',['jacobian_foot_left'],'vars',{z p});

matlabFunction(rE_r,'file',['position_foot_right'],'vars',{z p});
matlabFunction(drE_r,'file',['velocity_foot_right'],'vars',{z p});
%matlabFunction(ddrE_r,'file',['acceleration_foot_right'],'vars',{z p});
matlabFunction(J_r ,'file',['jacobian_foot_right'],'vars',{z p});


matlabFunction(dJ_l ,'file',['jacobian_dot_foot_left'],'vars',{z p});
matlabFunction(dJ_r ,'file',['jacobian_dot_foot_right'],'vars',{z p});

matlabFunction(Grav_Joint_Sp ,'file', ['Grav_walker'] ,'vars',{z p});
matlabFunction(Corr_Joint_Sp ,'file', ['Corr_walker']     ,'vars',{z p});
matlabFunction(F2Q(Fy*ihat(1:2), rE_l),'file',['Q_f_l'],'vars',{z p Fy});
matlabFunction(F2Q(Fy*ihat(1:2), rE_r),'file',['Q_f_r'],'vars',{z p Fy});

matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});

