%% stack_simulation.m
%C/A code transmission and receiver position acquisition simulation
%This function should serve as an example of how to use some of the
%toolkits funtions. 
%
% GPS STACK TOOLBOX
% Javier Antoran & Alberto Mur
% April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all
clc
% set env variables
ROOTDIR = fileparts(get_lib_path);
almFile = strcat(ROOTDIR,'/files/almanac/W918.alm');
ephFile = strcat(ROOTDIR,'/files/ephemeris/brdc0920.17n');

% read rinex and set contextual parameters
[eph, head] = read_rinex_nav(ephFile, 1:32);
[~, gps_sec] = cal2gpstime([2017 04 04 16 51 30]);
time = gps_sec + head.leapSeconds;
satp = eph2ecef(eph, time);

% Ellipsoid parameters
WGS84.a = 6378137;
WGS84.e2 = (8.1819190842622e-2).^2;

%define receiver position
Rpos = [ 3.894192036606761e+06 3.189618244369670e+05 5.024275884645306e+06]; % ECEF

%get visible SV
rcv_lla = [ 0 0 0 ];
[rcv_lla(1), rcv_lla(2), rcv_lla(3)] = xyz2lla(Rpos(1), Rpos(2), Rpos(3), WGS84.a, WGS84.e2);
e_mask = 30; %smallest angle between receiver and satellites for these to be visible 
visible_SV = visible_sv( satp, rcv_lla, e_mask );
%visible_SV = 1:32
[eph, head] = read_rinex_nav(ephFile, visible_SV);
%define constants
c = 2.99792458e8;
base_clock = 10.23e6;

%set sample frequency
L = 50; %samples per CA bit. fm = fchip * L. More renults ib better precision

% GPS signal generation

% generate C/A signals
sCA = CA_gen(L, visible_SV);

%generate doppler shift (use doppler function or create radom shift)
%f_vec = 2000 * randn(1, length(visible_SV)); % std is set to 2000 hz
%for better performance we can create our shift from a predefined set of values
shift_vals = (-1e4:5000:1e4);
f_vec = randi(length(shift_vals), 1, length(visible_SV));
doppler_carrier = generate_doppler( shift_vals(f_vec), L );
sCA_dop = sCA .* doppler_carrier; %add doppler shift

% pass signal through channel: distance + iono + tropo + clock drift +
% relativistic
[delay_CA, cicles, prop_delay_0, sat_clock_offset_0, sat_clock_rel_0, iono_T_0, trop_T_equiv_0] = gps_channel(head, eph, time, Rpos, sCA_dop, L);
srx = sum(delay_CA, 1);
%% receiver: GPS signal decoding

%set constants
f_chip = base_clock / 10;
fm = f_chip * L; %More results in better precision
Tm = 1/fm;
Lchip = 2 .^ 10 - 1;
Tchip = Lchip/f_chip;
samples_chip = Lchip * L;

%obtain SV postions:
% phase and time: adquisition

[ aquired, pr_delay_abs_samples, phase_delay ] = SV_CA_doppler_search( srx, L, 1e4, 5000);
%plot_CA_fi_search(  srx, 23, L, 1e4, 500 )

%%
ROOTDIR = fileparts(get_lib_path);
almFile = strcat(ROOTDIR,'/files/almanac/W918.alm');
ephFile = strcat(ROOTDIR,'/files/ephemeris/brdc0920.17n');

[eph, head] = read_rinex_nav(ephFile, aquired);

satp = eph2ecef(eph, time);
SVx = satp(2,:);
SVy = satp(3,:);
SVz = satp(4,:);
pSV = [SVx; SVy; SVz]';
%
pr_delay_samples = mod(pr_delay_abs_samples, samples_chip);
pr_delay = pr_delay_samples * Tm;

pseudo_range0 = (pr_delay +  cicles * Tchip)' * c;

%We will apply the following formula:
%pr0 = real_range + clock_drift + relativistic + iono + tropo

%variables for drift (iterative)
af = [eph.af0; eph.af1; eph.af2]';%bias, drift, drift rate
Toc = eph.toc; %time of clock of specific satellite relative to deltaT
% TOC IS UPDATED BY NAV MESSAGE i believe



% SAT CLOCK RELATIVISTIC TIME CORRECTION non iterative
A = eph.sqrtA.^2;
meanAnomaly  = eph.M0;
GM = 3.986005e14; %Grav constant * Earth mass
F = -2 * sqrt(GM) / c.^2;
tgd = eph.TGD; %group delay
T_amb=20; %degrees celsius
P_amb=101; %kilo pascal
P_vap=.86;
%delay because time in space is not the same as on earth
clock_relativistic = Error_Satellite_Clock_Relavastic(F,eph.e,A,meanAnomaly,tgd); %sec 
%%
R_rel_offset = clock_relativistic * c; 

pseudo_range0 = (pseudo_range0 - R_rel_offset);

% iterative decoding
rec_pos = [0 0 0];
%[G0, delta_x, N_rec_pos(1, :) ,B0]=Gen_G_DX_XYZ_B(pSV, rec_pos, pseudo_range0);
i = 1;

pseudo_range(1,:) = pseudo_range0;
%rec_pos = N_rec_pos(1, :);
delta_x = ones(1,4) * 1000;

while (norm(delta_x(1:3)) > 1)
   
    %clock drift
    Ttr = time - pseudo_range(i,:) ./ c;%gps time seconds
    clock_drift=Error_Satellite_Clock_Offset(af,Ttr,Toc); %sec
    R_c_offset(i, :) = clock_drift * c;
    
    %ionospheric delay
    Alpha = [head.A0 head.A1 head.A2 head.A3];
    Beta = [head.B0 head.B1 head.B2 head.B3];
    iono_T = Error_Ionospheric_Klobuchar(rec_pos, pSV ,Alpha, Beta, time);%(Sec)
    R_iono(i, :)=iono_T * c;   
    
    %tropospheric
    
    R_trop(i, :) = Error_Tropospheric_Hopfield(T_amb,P_amb,P_vap, rec_pos, pSV);
    
    % real range aprox
    
    pseudo_range(i + 1, :) = pseudo_range0 -  R_c_offset(i, :) - R_iono(i, :) - R_trop(i, :);
    
    i = i + 1;
    [G_mat, delta_x, N_rec_pos(i, :)]=iterate_pr2xyz(pSV, rec_pos, pseudo_range(i, :));
    rec_pos = N_rec_pos(i, :);
    clock_bias(i) = delta_x(4);
end

text = sprintf('Error, x: %d y: %d z: %d', (Rpos - rec_pos))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot position convergence
figure
subplot(3,1,1)
%hold on
plot(N_rec_pos(:,1))
line(xlim, [rec_pos(1),rec_pos(1)],'Color','r');
title('acquisition x')
xlabel('iterations')
ylabel('(m)')
legend('iteration value', 'final value')
grid on
subplot(3,1,2)
%hold on
plot(N_rec_pos(:,2))
title('acquisition y')
line(xlim, [rec_pos(2),rec_pos(2)],'Color','r');
xlabel('iterations')
ylabel('(m)')
legend('iteration value', 'final value')
grid on
subplot(3,1,3)
%hold on
plot(N_rec_pos(:,3))
title('acquisition z')
line(xlim, [rec_pos(3),rec_pos(3)],'Color','r');
xlabel('iterations')
ylabel('(m)')
legend('iteration value', 'final value')
grid on
%% pseudorange correction convergence
figure
index = aquired;
for i = 1:length(index)
    subplot(length(index), 1, i)
    plot(pseudo_range(:, i))
    line(xlim, [c*prop_delay_0(i), c*prop_delay_0(i)], 'Color', 'r');
    t = sprintf('pseudorange SV %d', index(i));
    title(t)
    xlabel('iterations')
    ylabel('pseudorange (m)')
    legend('calculated value', 'real value')
    grid on
end

%% pseudorrange error convergence **sat clock offset seems large**
figure()
N_SV = aquired(1);
subplot(3, 1, 1)
plot(R_c_offset(:, N_SV))
line(xlim, [c*sat_clock_offset_0(N_SV), c*sat_clock_offset_0(N_SV)], 'Color', 'r');
t = sprintf('clock drift SV %d', N_SV);
title(t)
xlabel('iterations')
ylabel('error (m)')
legend('calculated value', 'real value')
grid on
subplot(3, 1, 2)
plot(R_iono(:, N_SV))
line(xlim, [c*iono_T_0(N_SV), c*iono_T_0(N_SV)], 'Color', 'r');
t = sprintf('Ionospheric delay SV %d', N_SV);
title(t)
xlabel('iterations')
ylabel('error (m)')
legend('calculated value', 'real value')
grid on
subplot(3, 1, 3)
plot(R_trop(:, N_SV))
line(xlim, [c*trop_T_equiv_0(N_SV), c*trop_T_equiv_0(N_SV)], 'Color', 'r');
t = sprintf('tropospheric delay SV %d', N_SV);
title(t)
xlabel('iterations')
ylabel('error (m)')
legend('calculated value', 'real value')
grid on
%% Clock bias **B = deltax(4)**

%% G matrix and DOPs

Q = inv(G_mat' * G_mat);
GDOP=sqrt(trace(Q))
PDOP=sqrt(Q(1,1)+Q(2,2)+Q(3,3))
HDOP=sqrt(Q(1,1)+Q(2,2))
VDOP=sqrt(Q(1,1))
TDOP=sqrt(Q(4,4))