%%
clear
q0 = [0 0 0 0 0 0];
qf = [0 -pi/6 -pi/6 0 -pi/6 0];

t = linspace(0,3,50);
path_1 = zeros(3,50);
path_2 = quintic_trajectory(q0(2),0,0,qf(2),0,0,t);
path_3 = quintic_trajectory(q0(3),0,0,qf(3),0,0,t);
path_4 = zeros(3,50);
path_5 = quintic_trajectory(q0(5),0,0,qf(5),0,0,t);
path_6 = zeros(3,50);

Tau = generate_torque(path_1, path_2, path_3, path_4, path_5, path_6, 0);
figure;
tau_plot(t, Tau);

function tau_subplot(m, n, index, x, y)
ax = subplot(m,n,index);
plot(ax,x,y);
xlabel('Time (s)');
ylabel('Torque (N*m)');
title(['Joint #' num2str(index) ' Torque Output']);
end

function tau_plot(t,Tau)
sgtitle('Torque Controller Output for the 2-Link Robot');
max = size(Tau,1);
for p=1:max
    tau_subplot(3,2,p,t,Tau(p,:));
end
end

function plot_joints(path_1, path_2, t, title)
ax1 = subplot(3,2,1);
plot(ax1,t, path_1(1,:));
xlabel('Time (s)');
ylabel('Position (mm)');
ax2 = subplot(3,2,3);
plot(ax2,t, path_1(2,:));
xlabel('Time (s)');
ylabel('Velocity (mm/s)');
ax3 = subplot(3,2,5);
plot(ax3,t, path_1(3,:));
xlabel('Time (s)');
ylabel('Acceleration (mm/s^{2})');

ax1 = subplot(3,2,2);
plot(ax1,t, path_2(1,:));
xlabel('Time (s)');
ylabel('Position (mm)');
ax2 = subplot(3,2,4);
plot(ax2,t, path_2(2,:));
xlabel('Time (s)');
ylabel('Velocity (mm/s)');
ax3 = subplot(3,2,6);
plot(ax3,t, path_2(3,:));
xlabel('Time (s)');
ylabel('Acceleration (mm/s^{2})');

sgtitle(title);
end
