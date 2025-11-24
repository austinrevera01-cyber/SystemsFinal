function results = Systemsfinalprod(user_options)
clc;close all;
% If the caller provided an options struct, use it. Otherwise start fresh.
if nargin > 0
    opts = user_options;
else
    opts = struct();
end

function params = default_parameters()
    params = struct();
    params.encoder.counts_per_rev = 500 * 4;
    params.max_encoder  = 4096;
    params.motor.R      = 0.9;
    params.motor.L      = 0.45e-3;
    params.gear.N       = 21*15;
    params.vehicle.mass = 720;
    params.vehicle.Iz   = 1200;
    params.vehicle.a    = 1.72;
    params.vehicle.b    = 1.25;
    params.vehicle.Cf   = 100000;
    params.vehicle.Cr   = 100000;
    params.vehicle.Kt   = 0.0259; 
end

opts.Vel                = get_option(opts, 'Vel', 8);
opts.Ts                 = get_option(opts, 'Ts', 0.001);
opts.id_duration        = get_option(opts, 'id_duration', 3.5);
opts.id_voltage         = get_option(opts, 'id_voltage', 1);
opts.test_duration      = get_option(opts, 'test_duration', 3);

params = default_parameters();
if isfield(opts, 'params')
    params = merge_structs(params, opts.params);
end

%% Identification: collect step response from p-code to solve for Je and Be

clear pcode_identification;
id_time = (0:opts.Ts:opts.id_duration)'; % Time vector for identification
const_volts=zeros(size(id_time));
const_volts(:) = opts.id_voltage;
SS_values = pcode_identification(const_volts, id_time, opts, params);
const_rack_model = steering(params, SS_values);
%% Build the yaw-rate bicycle model and cascade with rack model
[yaw_tf] = build_bicycle_model(params, opts.Vel);
const_voltage_to_yaw_model = series(const_rack_model, yaw_tf);


figure('Name','Bode');
bode(const_voltage_to_yaw_model);
figure('Name','Eigenvalues');
p = pole(const_voltage_to_yaw_model);
plot(p, 0 , 'X', 'MarkerSize',10, 'LineWidth',3);
xlabel('Real'); ylabel('Imag');
title('Poles / Eigenvalues'); grid on;

figure('Name', 'Yaw Rate Comparison Vonstant Input');
plot(id_time, SS_values.yaw_rate, 'r', 'DisplayName', 'Yaw Rate'); hold on;
step(const_voltage_to_yaw_model,opts.id_duration);
grid on; xlabel('Time [s]'); ylabel('Yaw Rate [Rad/s]');
legend('show'); title('Yaw Rate Comparison');


%Run pcode with varriable V
clear collect_pcode_response;
id_time = (0:opts.Ts:opts.test_duration)'; % Time vector for identification
varr_volts = sin(5*id_time);
varr_data = collect_pcode_response(varr_volts, id_time, opts, params,SS_values);

figure('Name', 'Yaw Rate Comparison Varrible Input');
subplot(2,1,1)
plot(id_time,varr_volts, 'b', 'DisplayName', 'Input Voltage');
grid on; xlabel('Time [s]'); ylabel('Volts [V]');
subplot(2,1,2)
plot(id_time, varr_data.yaw_rate, 'r', 'DisplayName', 'P-code'); hold on;
plot(id_time, varr_data.sim_yaw, 'b', 'DisplayName', 'Transfer Function')
grid on; xlabel('Time [s]'); ylabel('Yaw Rate [Rad/s]');
legend('show'); title('Yaw Rate Comparison');

end


%%%%%%%%%%%%%%%%%%%
%%% Local Functions
%%%%%%%%%%%%%%%%%%%

function merged = merge_structs(base, override)
    merged = base;
    names = fieldnames(override);
    for k = 1:numel(names)
        f = names{k};
        if isstruct(override.(f))
            merged.(f) = merge_structs(base.(f), override.(f));
        else
            merged.(f) = override.(f);
        end
    end
end

function val = get_option(opts, name, default)
    if isfield(opts, name)
        val = opts.(name);
    else
        val = default;
    end
end