function controller = controller_dev(params,velocity,SS_values)

    N = params.gear.N;
    m = params.vehicle.mass;
    Iz = params.vehicle.Iz;
    a = params.vehicle.a;
    b = params.vehicle.b;
    Cf = params.vehicle.Cf;
    Cr = params.vehicle.Cr;
    Kt = params.vehicle.Kt;
    C0 = Cf + Cr;
    C1 = a*Cf + b*Cr;
    C2 = (a^2)*Cf + (b^2)*Cr;
    K = SS_values.K
    Je= SS_values.Je;
    Be= SS_values.Be;
    tau = Je/Be
    zeta=.7;
    TTS = .6;
    omega_n = 4 / (zeta*TTS)



    A = a*Cf / Iz
    B = (a+b)*Cr*Cf / (Iz*m*velocity)
    C = (C0*Iz + C2*m)/(Iz*m*velocity)
    D = ((C0*C2 - C1*m*velocity^2) - C1^2)/(Iz*m*velocity^2)

    controller.Kp1 = (2*tau*omega_n - 1)/K;
    controller.Ki1 = (tau*omega_n^2)/(K);

    omega_n2 = 4 / (zeta*(TTS/2))

    controller.Kp2 = (-A*C*(omega_n2^2) + 2*A*D*omega_n2*zeta - B*D + B*(omega_n2^2))/((A^2)*(omega_n2^2) - 2*A*B*omega_n2*zeta + B^2);
    controller.Kd2 = (A*D - A*omega_n2^2 - B*C + 2*B*omega_n2*zeta)/((A^2)*(omega_n2^2) - 2*A*B*omega_n2*zeta + B^2);
end