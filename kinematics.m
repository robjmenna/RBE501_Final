syms th1 th2 th3 th4 th5 th6 pi

% Define the DH parameters.
a = [50 330 35 0 0 0];
alpha = [pi/2 0 pi/2 -pi/2 pi/2 0];
d = [320 0 0 335 0 80];
theta = [th1 (th2+pi/2) th3 th4 th5 th6];

% Define the homogenous transformation matrix.
T10 = DH(alpha(1), a(1), d(1), theta(1));
T21 = DH(alpha(2), a(2), d(2), theta(2));
T32 = DH(alpha(3), a(3), d(3), theta(3));
T43 = DH(alpha(4), a(4), d(4), theta(4));
T54 = DH(alpha(5), a(5), d(5), theta(5));
Tf5 = DH(alpha(6), a(6), d(6), theta(6));

T20 = simplify(T10*T21);
T30 = simplify(T20*T32);
T40 = simplify(T30*T43);
T50 = simplify(T40*T54);
Tf0 = simplify(T50*Tf5);
display(Tf0);

% Plot the position of the robot.
theta_pos = [0 0 0 0 0 0];
T10_temp = subs(T10, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);
T20_temp = subs(T20, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T30_temp = subs(T30, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T40_temp = subs(T40, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T50_temp = subs(T50, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);
T60_temp = subs(Tf0, [th1 th2 th3 th4 th5 th6], [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]); 

x = [0 T10_temp(1,4) T20_temp(1,4) T30_temp(1,4) T40_temp(1,4) T50_temp(1,4) T60_temp(1,4)];
y = [0 T10_temp(2,4) T20_temp(2,4) T30_temp(2,4) T40_temp(2,4) T50_temp(2,4) T60_temp(2,4)];
z = [0 T10_temp(3,4) T20_temp(3,4) T30_temp(3,4) T40_temp(3,4) T50_temp(3,4) T60_temp(3,4)];

plot3(x,y,z,'-ro','LineWidth',4,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);

% Derive the Jacobian.
Jv = [diff(Tf0(1:3,4), th1) diff(Tf0(1:3,4), th2) diff(Tf0(1:3,4), th3) diff(Tf0(1:3,4), th4) diff(Tf0(1:3,4), th5) diff(Tf0(1:3,4), th6)];
Jv = simplify(Jv);
Jw = [get_z(T10) get_z(T20) get_z(T30) get_z(T40) get_z(T50) get_z(Tf0)];
Jw = simplify(Jw);
J = [Jv; Jw];
display(J);

% Return a transformation matrix based on the DH parameters.
function T = DH(alpha, a, d, theta)
T = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha)  a*cos(theta);
     sin(theta) cos(theta)*cos(alpha)  -cos(theta)*sin(alpha) a*sin(theta);
     0          sin(alpha)             cos(alpha)             d;
     0          0                      0                      1];
end

% Calculate an element of the angular component of the Jacobian.
function z = get_z(T)
    R = T(1:3,1:3);
    z = R * [0;0;1];
end