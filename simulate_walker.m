function simulate_walker()
    %% Definte fixed paramters
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.8; 

    ramp = 11*(pi/180); % in radians, (-) is uphill, (+) is downhill
    
    rest_coeff = 0.8;
    fr_coeff = 0.6; %0.7
    ground_height = -0.2;

    %% Parameter vector
    p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g ramp]'; 
    Fn = (m4)*g;
    K_c = 5000;
    D_c = 20;
    %% Simulation Parameters Set 2 -- Operational Space Control
    p_traj.omega_l = -10;
    p_traj.x_0_l   = 0;
    p_traj.y_0_l   = -0.2; %-.2;
    p_traj.r_l     = 0.02; %0.025;
    p_traj.offset = pi;

    p_traj.omega_r = -10;
    p_traj.x_0_r   = 0; %-.02;
    p_traj.y_0_r   = -0.2; %-.12;
    p_traj.r_r     = 0.02;

    % for ellipsoid stuff
    p_traj.a = 2;
    p_traj.b = 1;
    p_traj.theta = 0.2;%0.1;
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 10;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [0;0;pi/5; pi/4; pi/6; pi/4; 0;0;0;0; 0;0];
    z_out = zeros(12,num_step);
    z_out(:,1) = z0;
    
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p, p_traj, fr_coeff, Fn, rest_coeff, K_c, D_c);
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        z_out(1:6,i+1) = z_out(1:6,i) + z_out(7:end,i+1)*dt;
    end

    %% Animate Solution
    figure(6); clf;
    hold on
    
    % Ground Q2.3
    plot([-.2 4],[ground_height-.2*sin(-ramp) ground_height+4*sin(-ramp)],'k'); 
    
    animateSol(tspan, z_out,p, p_traj, ramp);

end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 1500; % Spring stiffness X
    K_y = 1500; % Spring stiffness Y
    D_x = 10;  % Damping X
    D_y = 10;  % Damping Y

    % Desired position of foot is a circle
    omega_swing_l = p_traj.omega_l;
    r0 = r0_walker(z,p);
    r0 = r0(1:2);
    dr0 = dr0_walker(z,p);
    dr0 = dr0(1:2);

    a = p_traj.a; 
    b = p_traj.b;
    theta = p_traj.theta;
    rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rEd_l = r0 +  [p_traj.x_0_l; p_traj.y_0_l] + ...
        rotation_matrix*p_traj.r_l * [a*cos(omega_swing_l * t + p_traj.offset); b*sin(omega_swing_l * t + p_traj.offset)];

    vEd_l = dr0 + rotation_matrix*p_traj.r_l*[-a*sin(omega_swing_l*t+ p_traj.offset)*omega_swing_l    ...
                     b*cos(omega_swing_l*t+ p_traj.offset)*omega_swing_l ]';
    % Actual position and velocity 
    rE_l = position_foot_left(z,p);
    vE_l = velocity_foot_left(z,p);

    omega_swing_r = p_traj.omega_r;

    rEd_r = r0 + [p_traj.x_0_r; p_traj.y_0_r] + ...
        rotation_matrix*p_traj.r_r * [a*cos(omega_swing_r * t); b*sin(omega_swing_r * t)];


    % Compute desired velocity & acceleration of foot
    vEd_r = dr0 + rotation_matrix*p_traj.r_r*[-a*sin(omega_swing_r*t)*omega_swing_r    ...
                     b*cos(omega_swing_r*t)*omega_swing_r ]';
  
    
    % Actual position and velocity 
    rE_r = position_foot_right(z,p);
    vE_r = velocity_foot_right(z,p);

    %%% Prev controller to op. space
    
    pVprEl = [K_x*(rE_l(1)-rEd_l(1)); K_y*(rE_l(2)- rEd_l(2))];
    vel_l = velocity_foot_left(z, p);
    vels = vel_l(1:2,:);
    F_dampl = [-D_x*(vels(1,:) - vEd_l(1)); -D_y*(vels(2,:) - vEd_l(2))];
    J_l  = jacobian_foot_left(z,p);
    taut = -transpose(J_l(1:2,:)) * pVprEl + transpose(J_l(1:2,:))*F_dampl;

    pVprEr = [K_x*(rE_r(1)-rEd_r(1)); K_y*(rE_r(2)- rEd_r(2))];
    vel_r = velocity_foot_right(z, p);
    vels = vel_r(1:2,:);
    F_dampr = [-D_x*(vels(1,:) - vEd_r(1)); -D_y*(vels(2,:) - vEd_r(2))];
    J_r  = jacobian_foot_right(z,p);
    taut2 = -transpose(J_r(1:2,:)) * pVprEr + transpose(J_r(1:2,:))*F_dampr;

    taut_both = taut + taut2;

    tau = taut_both(3:end);

end


function dz = dynamics(t,z,p,p_traj, fr_coeff, Fn, rest_coeff, K_c, D_c)
    % Get mass matrix
    A = A_walker(z,p);
    
    % Compute Controls
   tau = control_law(t,z,p,p_traj);
   %tau = [0;0;0;0];
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_walker(z,tau,p);

    [Fc_left, Fc_right] = contact_force(z, p, K_c, D_c, -0.2, fr_coeff, Fn, rest_coeff);

    % Compute the contribution of the contact force to the generalied force
    J_r = jacobian_foot_right(z,p);
    J_l = jacobian_foot_left(z,p);

    QF_R = J_r'*(Fc_right);    % contact force contributions to generalized forces
    QF_L = J_l'*(Fc_left);


%     QF_L = zeros(6, 1); % for testing
%     QF_R = zeros(6, 1); % for testing
    
    % Solve for qdd.
    qdd = A\(b + QF_L + QF_R);

    %qdd = A\(b);
    dz = 0*z;

    % Form dz
    dz(1:6) = z(7:12);
    dz(7:12) = qdd;
end

function [Fc_left, Fc_right] = contact_force(z, p, K_c, D_c, yC, fr_coeff, Fn, rest_coeff)
% left leg
    dq_l = z(7:end);
    rE_l = position_foot_left(z,p);
    drE_l = velocity_foot_left(z,p);
    J_l = jacobian_foot_left(z,p);
    M = A_walker(z,p);
    C_l = rE_l(2) - yC;
    dC_l = drE_l(2);
    Fc_left = zeros(2,1);
    Fc_right = zeros(2,1);
   
    if C_l > 0 || (-K_c*C_l  - D_c*dC_l) < 0
        Fc_left(2) = 0;
        Fc_left(1) = 0;
    else
        Fc_left(2) = (-K_c*C_l  - D_c*dC_l);
        Fc_left(1) = Fc_left(2)*fr_coeff*sign(-drE_l(1)); 
    end
% right leg 
    dq_r = z(7:end);
    rE_r = position_foot_right(z,p);
    drE_r = velocity_foot_right(z,p);
    C_r = rE_r(2) - yC;
    dC_r = drE_r(2);
    if C_r > 0 || (-K_c*C_r  - D_c*dC_r) < 0
        Fc_right(2) = 0;
        Fc_right(1) = 0;
    else
        Fc_right(2) = (-K_c*C_r  - D_c*dC_r); 
        Fc_right(1) = Fc_right(2)*fr_coeff*sign(-drE_r(1));
    end

end

function animateSol(tspan, x,p, p_traj, ramp)
%     x = z_out
%     myVideo = VideoWriter('myVideoFile'); %open video file
%     myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
%     open(myVideo)

%     left_foot_traj = plot([0],[0],'LineWidth',1);
%     right_foot_traj = plot([0],[0],'LineWidth',1);

    h_OB_l = plot([0],[0],'LineWidth',2);
    h_AC_l = plot([0],[0],'LineWidth',2);
    h_BD_l = plot([0],[0],'LineWidth',2);
    h_CE_l = plot([0],[0],'LineWidth',2);

    h_OB_r = plot([0],[0],'LineWidth',2);
    h_AC_r = plot([0],[0],'LineWidth',2);
    h_BD_r = plot([0],[0],'LineWidth',2);
    h_CE_r = plot([0],[0],'LineWidth',2);

    R = [cos(-ramp), -sin(-ramp); sin(-ramp), cos(-ramp)];
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    %axis([-.2 .2 -.3 .1]);
    axis([-.1 4 -1 .1]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_walker(z,p);

        r0 = keypoints(:,1); % x and y of body
        rA_l = R*keypoints(:,2); 
        rB_l = R*keypoints(:,3);
        rC_l = R*keypoints(:,4); 
        rD_l = R*keypoints(:,5);
        rE_l = R*keypoints(:,6);

        % Desired left foot trajectory
        omega_swing_l = p_traj.omega_l;
        omega_swing_r = p_traj.omega_r;

        a = p_traj.a; 
        b = p_traj.b;
        theta = p_traj.theta;

        rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

        rEd_l = r0 +  [p_traj.x_0_l; p_traj.y_0_l] + ...
            rotation_matrix*p_traj.r_l * [a*cos(omega_swing_l * t + p_traj.offset); b*sin(omega_swing_l * t + p_traj.offset)];

        rEd_r = r0 + [p_traj.x_0_r; p_traj.y_0_r] + ...
            rotation_matrix*p_traj.r_r * [a*cos(omega_swing_r * t); b*sin(omega_swing_r * t)];

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title

        
        rEd_l = R*rEd_l;
        rEd_r = R*rEd_r;

        r0 = R*keypoints(:,1); %now rotate body coordinates
        
        % Plot desired left foot trajectory
        plot(rEd_l(1), rEd_l(2), 'ro', 'MarkerSize',1);
        plot(rEd_r(1), rEd_r(2), 'bo','MarkerSize',1);

        set(h_OB_l,'XData',[r0(1) rB_l(1)]);
        set(h_OB_l,'YData',[r0(2) rB_l(2)]);
        
        set(h_AC_l,'XData',[rA_l(1) rC_l(1)]);
        set(h_AC_l,'YData',[rA_l(2) rC_l(2)]);
        
        set(h_BD_l,'XData',[rB_l(1) rD_l(1)]);
        set(h_BD_l,'YData',[rB_l(2) rD_l(2)]);
        
        set(h_CE_l,'XData',[rC_l(1) rE_l(1)]);
        set(h_CE_l,'YData',[rC_l(2) rE_l(2)]);


        rA_r = R*keypoints(:,7);
        rB_r = R*keypoints(:,8);
        rC_r = R*keypoints(:,9); 
        rD_r = R*keypoints(:,10);
        rE_r = R*keypoints(:,11);
        
        set(h_OB_r,'XData',[r0(1) rB_r(1)]);
        set(h_OB_r,'YData',[r0(2) rB_r(2)]);
        
        set(h_AC_r,'XData',[rA_r(1) rC_r(1)]);
        set(h_AC_r,'YData',[rA_r(2) rC_r(2)]);
        
        set(h_BD_r,'XData',[rB_r(1) rD_r(1)]);
        set(h_BD_r,'YData',[rB_r(2) rD_r(2)]);
        
        set(h_CE_r,'XData',[rC_r(1) rE_r(1)]);
        set(h_CE_r,'YData',[rC_r(2) rE_r(2)]);

%         frame = getframe(gcf); %get frame
%         writeVideo(myVideo, frame);

        pause(0.01)
    end

%     close(myVideo)
end