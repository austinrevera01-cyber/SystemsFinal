function controller = controller_dev(params, velocity, SS_values, plot_opts)
%CONTROLLER_DEV Design and evaluate cascaded yaw control loops (Part B).
%   controller = CONTROLLER_DEV(params, velocity, SS_values, plot_opts)
%   designs PI/PD gains for steering angle, yaw rate, and heading loops and
%   exercises the closed-loop responses required for Part B of the project.
%
%   Inputs:
%       params      - vehicle parameter struct (mass, geometry, cornering)
%       velocity    - longitudinal velocity [m/s] used for linearization
%       SS_values   - struct containing rack gain (K), inertia (Je),
%                     viscous damping (Be) from identification
%       plot_opts   - optional struct with logical fields:
%                       .show_plots (default true)
%                       .figure_prefix (string used to group figures)
%
%   Outputs:
%       controller  - struct containing tuned gains, closed-loop transfer
%                     functions, stability metrics, and steady-state tests

    if nargin < 4
        plot_opts = struct();
    end
    if ~isfield(plot_opts, 'show_plots')
        plot_opts.show_plots = true;
    end
    if ~isfield(plot_opts, 'figure_prefix')
        plot_opts.figure_prefix = 'Part B - ';
    end

    % Vehicle parameters
    m  = params.vehicle.mass;
    Iz = params.vehicle.Iz;
    a  = params.vehicle.a;
    b  = params.vehicle.b;
    Cf = params.vehicle.Cf;
    Cr = params.vehicle.Cr;

    % Rack parameters from identification
    K   = SS_values.K;
    Je  = SS_values.Je;
    Be  = SS_values.Be;
    tau = Je / Be;

    % Desired dynamics
    zeta = 0.7;   % damping ratio
    TTS  = 0.6;   % target time-to-settle [s]

    % Bicycle model shorthand
    C0 = Cf + Cr;
    C1 = a*Cf + b*Cr;
    C2 = (a^2)*Cf + (b^2)*Cr;

    A = a*Cf / Iz;
    B = (a+b)*Cr*Cf / (Iz*m*velocity);
    C = (C0*Iz + C2*m)/(Iz*m*velocity);
    D = ((C0*C2 - C1*m*velocity^2) - C1^2)/(Iz*m*velocity^2);

    % Inner steering loop (PI)
    omega_n_delta = 4 / (zeta*TTS/6);
    controller.Kp1 = (2*tau*omega_n_delta - 1) / K;
    controller.Ki1 = (tau*omega_n_delta^2) / K;

    % Yaw-rate loop (PD)
    omega_n_r = 4 / (zeta*(TTS/3));
    controller.Kp2 = -((-A*C*(omega_n_r^2) + 2*A*D*omega_n_r*zeta - B*D + B*(omega_n_r^2)) ...
                       / ((A^2)*(omega_n_r^2) - 2*A*B*omega_n_r*zeta + B^2)) + 5;
    controller.Kd2 = -((A*D - A*omega_n_r^2 - B*C + 2*B*omega_n_r*zeta) ...
                       / ((A^2)*(omega_n_r^2) - 2*A*B*omega_n_r*zeta + B^2));

    % Heading loop (PI)
    omega_n_psi = 4 / (zeta*TTS);
    controller.Kp3 = omega_n_psi*2*zeta + 2;
    controller.Ki3 = omega_n_psi^2;

    s = tf('s');

    % Plants
    G_delta  = K/(tau*s + 1);              % V -> steering angle
    G_rdelta = (A*s + B)/(s^2 + C*s + D);  % steering -> yaw rate

    % Controllers
    C_delta = controller.Kp1 + controller.Ki1/s;    % inner PI (δ-loop)
    C_r     = controller.Kp2 + controller.Kd2*s;    % yaw-rate PD
    C_psi   = controller.Kp3 + controller.Ki3/s;    % outer PI (heading)

    %% Closed-loop interconnections
    % 1) Inner steering loop: δ_ref -> δ
    T_delta = feedback(C_delta*G_delta, 1);

    % 2) Yaw-rate loop: r_ref -> r (uses closed δ-loop as actuator)
    L_r = C_r * G_rdelta * T_delta;
    T_r = feedback(L_r, 1);

    % 3) Heading loop: ψ_ref -> ψ
    G_psi_eff = T_r / s;             % r_ref -> ψ is yaw-loop then integrator
    L_psi = C_psi * G_psi_eff;
    T_psi = feedback(L_psi, 1);

    %% Test metrics (Part B requirements)
    controller.loops.delta.tf = minreal(T_delta);
    controller.loops.delta.poles = pole(T_delta);
    controller.loops.delta.step = stepinfo(T_delta);
    controller.loops.delta.dc_gain = dcgain(T_delta);

    controller.loops.r.tf = minreal(T_r);
    controller.loops.r.poles = pole(T_r);
    controller.loops.r.step = stepinfo(T_r);
    controller.loops.r.bandwidth = bandwidth(T_r);

    controller.loops.psi.tf = minreal(T_psi);
    controller.loops.psi.poles = pole(T_psi);
    controller.loops.psi.step = stepinfo(T_psi);
    controller.loops.psi.dc_gain = dcgain(T_psi);

    controller.TF = controller.loops.psi.tf;

    if plot_opts.show_plots
        prefix = plot_opts.figure_prefix;
        quick_step_plot(T_delta, [prefix 'Inner steering: \delta_{ref} -> \delta']);
        quick_step_plot(T_r, [prefix 'Yaw-rate loop: r_{ref} -> r']);
        quick_step_plot(T_psi, [prefix 'Heading loop: \psi_{ref} -> \psi']);

        pole_plot(controller.loops.delta.poles, [prefix 'Steering loop poles']);
        pole_plot(controller.loops.r.poles, [prefix 'Yaw-rate loop poles']);
        pole_plot(controller.loops.psi.poles, [prefix 'Heading loop poles']);

        figure('Name', [prefix 'Heading loop Bode']);
        bode(T_psi);
        grid on;
    end
end

%%%%%%%%%%%%%%%%%%%
%%% Local helpers
%%%%%%%%%%%%%%%%%%%
function quick_step_plot(T, title_str)
    figure('Name', title_str);
    step(T);
    grid on;
    title(title_str);
end

function pole_plot(poles, title_str)
    figure('Name', title_str);
    plot(real(poles), imag(poles), 'x', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('Real'); ylabel('Imag'); grid on; title(title_str);
end