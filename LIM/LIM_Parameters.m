clear all;
clc;
%% Compute the Lorentz force required once known the POD's mass and velocity desired

% Funtion inputs
mass=200;                   % [kg]
max_vel= 120;               % [km/h]
track_length = 150;         % [m]

acc_perc= 0.6; % Percentage of the total distance used for accelerate
dec_perc= 0.4; % Percentage of the total distance used for decelerate

% Velocity, acceleration and deceleration

v_max = max_vel/3.6;                        % [m/s]

acc_dist = acc_perc*track_length;           % [m]
a_acc = v_max^2/(2*acc_dist);               % [m/s^2]
g_acc = a_acc/9.81;

dec_dist = dec_perc*track_length;
a_dec = -(v_max^2)/(2*dec_dist);
g_dec = a_dec/9.81;

% Aplying Newton's third law we obtain the thrust 

F_Lorentz_acc = mass*a_acc;                     % [N]
F_Lorentz_dec = mass*a_dec;                     % [N]

% Trajectory profile

run_time1 = v_max/a_acc;                                    % Acceleration
run_time2 = (1-acc_perc-dec_perc)*track_length/v_max;       % Max speed
run_time3 = -v_max/a_dec;                                   % Deceleration                                   

total_time = run_time1 + run_time2 + run_time3;             % [seconds]

t = 1;
h=0.05; % Time step [s]

while ((t)*h)<run_time1
    current_vel(t) = a_acc*(h*t);
    current_dist(t) = 0.5*a_acc*(h*t)^2;
    t=t+1;
end
while ((t)*h)<(run_time1+run_time2)
    current_vel(t) = v_max;
    current_dist(t) = acc_dist+v_max*((t*h)-run_time1);
    t=t+1;
end
while ((t)*h)<total_time
    current_vel(t) = v_max+a_dec*((h*t)-run_time1-run_time2);
    current_dist(t) = (track_length-dec_dist)+v_max*((h*t)-run_time1-run_time2)+0.5*a_dec*((h*t)-run_time1-run_time2)^2;
    t=t+1;
end

plot(current_dist(:),current_vel(:))

%% LIM sizing

max_width = 92.54;      % [mm]
number_poles = 2;
number_slots = 24;
spp = 4;                % Slots per phase 
CI = 140*sqrt(2);       % Peak current intensity per phase [A] 
turns = 25;
tpp = turns*spp;        % Turns per phase

railway_width = 10.46; % Given by the competition [mm]
distance2railway = 3.5; % Can be modified to increase the thrust [mm]
air_gap = railway_width+2*distance2railway;  % [mm]

% Performance constants
slip = 0.2;                       % (vs-vr)/vs
vr = v_max*1e3;                   % Mechanical rotor speed - Real speed [mm/s]
vs = vr/(1-slip);                 % Synchronous rotor speed - Ideal limit that wont be reached [mm/s]
cond_rotor = 2.5e4;               % Rotor conductivity of Aluminium 6061 T6 [S/mm]
d_rotor = railway_width;          % Rotor thickness [mm]
mu = 5000;                           % Magnetic permeability [H/mm]
MMF = turns * spp * CI;           % Magnetomotive Force [At]

% We will compute the thrust for diferent pole pitches [mm] 

for p_p = 1:500

    % Amplitude of the magnetic flux constant
    B1(p_p) = (3*sqrt(2)*MMF)/sqrt((pi*air_gap/mu)^2+(vs-vr)^2*(p_p*d_rotor*cond_rotor)^2);
    
    % Phase Shift
    phase_shift = atan((pi()*air_gap)/(cond_rotor*d_rotor*p_p*mu*(vs-vr)));
    fun = @(y,t) real(-cond_rotor*(vs-vr)*B1(p_p).*exp(1j*((pi/p_p)*(y-t*vs)+phase_shift))*B1(p_p).*exp(-1j*((pi/p_p)*(y-t*vs)+phase_shift)));

    Fint = integral(@(y) fun(y,100),0,2*p_p);

    Thrust(p_p) = max_width * d_rotor * Fint;

end

plot((1:500), Thrust(:))

% Something is wrong when computing the thrust as B1 seems to be
% wrong. We will continue with a pole pitch around 400mm. 

p_p = 400; 

% FREQUENCIES AND SKIN DEPTH

i=1;
for vel = 1:vr
    freq(i) = vel/(2*p_p); 
    skin_depth(i) = sqrt(1/(freq(i)))*100.49;  % [mm]
    i=i+1;
end

%plot(freq(:),skin_depth(:))


% The two stators operate at their maximum when acting together. This 
% requires the magnetic field to penetrate through the I-beam and reach 
% the other stator. This is fulfilled if at all operating frequencies the 
% skin depth is GREATER than the track width. 
