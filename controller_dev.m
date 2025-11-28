function controller = controller_dev(params,velocity,SS_values)

    m = params.vehicle.mass;
    Iz = params.vehicle.Iz;
    a = params.vehicle.a;
    b = params.vehicle.b;
    Cf = params.vehicle.Cf;
    Cr = params.vehicle.Cr;
    C0 = Cf + Cr;
    C1 = a*Cf + b*Cr;
    C2 = (a^2)*Cf + (b^2)*Cr;
    K = SS_values.K;
    Je= SS_values.Je;
    Be= SS_values.Be;
    tau = Je/Be;
    zeta=.7;
    TTS = .6;
    omega_n = 4 / (zeta*TTS/6);



    A = a*Cf / Iz;
    B = (a+b)*Cr*Cf / (Iz*m*velocity);
    C = (C0*Iz + C2*m)/(Iz*m*velocity);
    D = ((C0*C2 - C1*m*velocity^2) - C1^2)/(Iz*m*velocity^2);

    controller.Kp1 = (2*tau*omega_n - 1)/K;
    controller.Ki1 = (tau*omega_n^2)/(K);

    omega_n2 = 4 / (zeta*(TTS/3));

    controller.Kp2 = -((-A*C*(omega_n2^2) + 2*A*D*omega_n2*zeta - B*D + B*(omega_n2^2))/((A^2)*(omega_n2^2) - 2*A*B*omega_n2*zeta + B^2)) + 5;
    controller.Kd2 = -((A*D - A*omega_n2^2 - B*C + 2*B*omega_n2*zeta)/((A^2)*(omega_n2^2) - 2*A*B*omega_n2*zeta + B^2));

    omega_n3 = 4 / (zeta*TTS);

    controller.Kp3 = omega_n3*2*zeta  + 2;
    controller.Ki3 = omega_n3^2;

    s = tf('s');

    % Plants
    G_delta  = K/(tau*s + 1);              % V -> steering angle
    G_rdelta = (A*s + B)/(s^2 + C*s + D);  % steering -> yaw rate
    
    % Controllers (use your numeric gains)
    C_delta = controller.Kp1 + controller.Ki1/s;    % inner PI (δ-loop)
    C_r     = controller.Kp2 + controller.Kd2*s;    % yaw-rate PD  <<< NOTE: *s*, not 1/s
    C_psi   = controller.Kp3 + controller.Ki3/s;    % outer PI (heading loop)
    
    %% 1) Inner steering loop: δ_ref -> δ
    T_delta = feedback(C_delta*G_delta, 1);        % negative feedback
    
    %% 2) Yaw-rate loop: r_ref -> r (uses closed δ-loop as actuator)
    L_r = C_r * G_rdelta * T_delta;
    T_r = feedback(L_r, 1);
    
    %% 3) Heading loop: ψ_ref -> ψ
    G_psi_eff = T_r / s;             % r_ref -> ψ is yaw-loop then integrator
    L_psi = C_psi * G_psi_eff;
    T_psi = feedback(L_psi, 1);      % final closed-loop from ψ_ref to ψ

    %%Test Values on controller gains
    %a1 = 1 + A*controller.Kd2;
    %a2 = C + A*controller.Kp2 + B*controller.Kd2;
    %a3 = D + B*controller.Kp2;
    
    controller.TF = minreal(T_psi);
    figure; step(T_delta); title('Inner steering: \delta_{ref} -> \delta');
    figure; step(T_r);     title('Yaw-rate loop: r_{ref} -> r');
    figure; step(T_psi);   title('Heading loop: \psi_{ref} -> \psi');
    pole(T_delta)
    pole(T_r)
    pole(T_psi)


    

end