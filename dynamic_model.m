%%
% This script will define the forward kinematics, the Jacobian, and the
% dynamic model for the Fanuc 6 d.o.f. robot and print the results in terms
% of the declared symbolic variables.

clear;
% The position of each joint.
syms q1 q2 q3 q4 q5 q6 pi;

% The inertia of the last three links.
syms I4 I5 I6

% The velocity of each joint.
syms dq1 dq2 dq3 dq4 dq5 dq6

% The acceleration of each joint.
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6

% The mass of each link.
syms m1 m2 m3 m4 m5 m6 g

q = [q1 q2 q3 q4 q5 q6];
dq = [dq1 dq2 dq3 dq4 dq5 dq6];
ddq = [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6];

% Define the link lengths as specified in the datasheet.
a1 = 50;
a2 = 330;
a3 = 35;
d1 = 320;
d4 = 335;
d6 = 80;

% Define the DH parameters.
a = [a1 a2 a3 0 0 0];
alpha = [pi/2 0 pi/2 -pi/2 pi/2 0];
d = [d1 0 0 d4 0 d6];
theta = [q1 (q2+pi/2) q3 q4 q5 q6];

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
T10_temp = subs(T10, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);
T20_temp = subs(T20, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T30_temp = subs(T30, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T40_temp = subs(T40, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);   
T50_temp = subs(T50, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]);
T60_temp = subs(Tf0, q, [theta_pos(1) theta_pos(2) theta_pos(3) theta_pos(4) theta_pos(5) theta_pos(6)]); 

x = [0 T10_temp(1,4) T20_temp(1,4) T30_temp(1,4) T40_temp(1,4) T50_temp(1,4) T60_temp(1,4)];
y = [0 T10_temp(2,4) T20_temp(2,4) T30_temp(2,4) T40_temp(2,4) T50_temp(2,4) T60_temp(2,4)];
z = [0 T10_temp(3,4) T20_temp(3,4) T30_temp(3,4) T40_temp(3,4) T50_temp(3,4) T60_temp(3,4)];

plot3(x,y,z,'-ro','LineWidth',4,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10);

% Derive the Jacobian.
Jv1 = simplify(jacobian(T10(1:3,4), q1));
Jw1 = simplify(get_z(T10));

Jv2 = simplify(jacobian(T20(1:3,4),[q1 q2]));
Jw2 = simplify([get_z(T10) get_z(T20)]);

Jv3 = simplify(jacobian(T30(1:3,4),[q1 q2 q3]));
Jw3 = simplify([get_z(T10) get_z(T20) get_z(T30)]);

Jv4 = simplify(jacobian(T40(1:3,4),[q1 q2 q3 q4]));
Jw4 = simplify([get_z(T10) get_z(T20) get_z(T30) get_z(T40)]);

Jv5 = simplify(jacobian(T50(1:3,4),[q1 q2 q3 q4 q5]));
Jw5 = simplify([get_z(T10) get_z(T20) get_z(T30) get_z(T40) get_z(T50)]);

Jvf = simplify(jacobian(Tf0(1:3,4),[q1 q2 q3 q4 q5 q6]));
Jwf = simplify([get_z(T10) get_z(T20) get_z(T30) get_z(T40) get_z(T50) get_z(Tf0)]);
display(Jvf);
display(Jwf);

%Find the Lagrangian by calculating the kinetic energy of each link
K(1) = trans_kinetic(dq1, Jv1, m1);
K(2) = trans_kinetic([dq1 dq2], Jv2, m2);
K(3) = trans_kinetic([dq1 dq2 dq3], Jv3, m3);
K(4) = trans_kinetic([dq1 dq2 dq3 dq4], Jv4, m4);
K(5) = trans_kinetic([dq1 dq2 dq3 dq4 dq5], Jv5, m5);
K(6) = trans_kinetic([dq1 dq2 dq3 dq4 dq5 dq6], Jvf, m6);

K(7) = rot_kinetic([dq1 dq2 dq3 dq4], Jw4, T40(1:3,1:3), I4);
K(8) = rot_kinetic([dq1 dq2 dq3 dq4 dq5], Jw5, T50(1:3,1:3), I5);
K(9) = rot_kinetic([dq1 dq2 dq3 dq4 dq5 dq6], Jwf, Tf0(1:3,1:3), I6);

Kt = 0;
for i = K
    Kt = Kt+i;
end

% Find the potential energy of each link.
P(1) = T10(3,4)*m1*g;
P(2) = T20(3,4)*m2*g;
P(3) = T30(3,4)*m3*g;
P(4) = T40(3,4)*m4*g;
P(5) = T50(3,4)*m5*g;
P(6) = Tf0(3,4)*m6*g;

Pt = 0;
for i = P
    Pt = Pt + i;
end

% The Lagrange equation
L = simplify(Kt-Pt);

A1 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 1);
A2 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 2);
A3 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 3);
A4 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 4);
A5 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 5);
A6 = get_A(L, [q1 q2 q3 q4 q5 q6], [dq1 dq2 dq3 dq4 dq5 dq6], [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], 6);

B1 = diff(L, q1);
B2 = diff(L, q2);
B3 = diff(L, q3);
B4 = diff(L, q4);
B5 = diff(L, q5);
B6 = diff(L, q6);

% Find the dynamic model of each link
Tau_l1 = A1 - B1;
Tau_l2 = A2 - B2;
Tau_l3 = A3 - B3;
Tau_l4 = A4 - B4;
Tau_l5 = A5 - B5;
Tau_l6 = A6 - B6;

% Find the elements of the Inertia Matrix
M11 = simplify(Tau_l1 - subs(Tau_l1,ddq1,0)) / ddq1;
M12 = simplify(Tau_l1 - subs(Tau_l1,ddq2,0)) / ddq2;
M13 = simplify(Tau_l1 - subs(Tau_l1,ddq3,0)) / ddq3;
M14 = simplify(Tau_l1 - subs(Tau_l1,ddq4,0)) / ddq4;
M15 = simplify(Tau_l1 - subs(Tau_l1,ddq5,0)) / ddq5;
M16 = simplify(Tau_l1 - subs(Tau_l1,ddq6,0)) / ddq6;

M21 = simplify(Tau_l2 - subs(Tau_l2,ddq1,0)) / ddq1;
M22 = simplify(Tau_l2 - subs(Tau_l2,ddq2,0)) / ddq2;
M23 = simplify(Tau_l2 - subs(Tau_l2,ddq3,0)) / ddq3;
M24 = simplify(Tau_l2 - subs(Tau_l2,ddq4,0)) / ddq4;
M25 = simplify(Tau_l2 - subs(Tau_l2,ddq5,0)) / ddq5;
M26 = simplify(Tau_l2 - subs(Tau_l2,ddq6,0)) / ddq6;

M31 = simplify(Tau_l3 - subs(Tau_l3,ddq1,0)) / ddq1;
M32 = simplify(Tau_l3 - subs(Tau_l3,ddq2,0)) / ddq2;
M33 = simplify(Tau_l3 - subs(Tau_l3,ddq3,0)) / ddq3;
M34 = simplify(Tau_l3 - subs(Tau_l3,ddq4,0)) / ddq4;
M35 = simplify(Tau_l3 - subs(Tau_l3,ddq5,0)) / ddq5;
M36 = simplify(Tau_l3 - subs(Tau_l3,ddq6,0)) / ddq6;

M41 = simplify(Tau_l4 - subs(Tau_l4,ddq1,0)) / ddq1;
M42 = simplify(Tau_l4 - subs(Tau_l4,ddq2,0)) / ddq2;
M43 = simplify(Tau_l4 - subs(Tau_l4,ddq3,0)) / ddq3;
M44 = simplify(Tau_l4 - subs(Tau_l4,ddq4,0)) / ddq4;
M45 = simplify(Tau_l4 - subs(Tau_l4,ddq5,0)) / ddq5;
M46 = simplify(Tau_l4 - subs(Tau_l4,ddq6,0)) / ddq6;

M51 = simplify(Tau_l5 - subs(Tau_l5,ddq1,0)) / ddq1;
M52 = simplify(Tau_l5 - subs(Tau_l5,ddq2,0)) / ddq2;
M53 = simplify(Tau_l5 - subs(Tau_l5,ddq3,0)) / ddq3;
M54 = simplify(Tau_l5 - subs(Tau_l5,ddq4,0)) / ddq4;
M55 = simplify(Tau_l5 - subs(Tau_l5,ddq5,0)) / ddq5;
M56 = simplify(Tau_l5 - subs(Tau_l5,ddq6,0)) / ddq6;

M61 = simplify(Tau_l6 - subs(Tau_l6,ddq1,0)) / ddq1;
M62 = simplify(Tau_l6 - subs(Tau_l6,ddq2,0)) / ddq2;
M63 = simplify(Tau_l6 - subs(Tau_l6,ddq3,0)) / ddq3;
M64 = simplify(Tau_l6 - subs(Tau_l6,ddq4,0)) / ddq4;
M65 = simplify(Tau_l6 - subs(Tau_l6,ddq5,0)) / ddq5;
M66 = simplify(Tau_l6 - subs(Tau_l6,ddq6,0)) / ddq6;

% Find the elements of the gravity matrix.
G11 = subs(Tau_l1, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);
G21 = subs(Tau_l2, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);
G31 = subs(Tau_l3, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);
G41 = subs(Tau_l4, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);
G51 = subs(Tau_l5, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);
G61 = subs(Tau_l6, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6],...
    [0 0 0 0 0 0 0 0 0 0 0 0]);

M = [M11 M12 M13 M14 M15 M16;
     M21 M22 M23 M24 M25 M26;
     M31 M32 M33 M34 M35 M36;
     M41 M42 M43 M44 M45 M46;
     M51 M52 M53 M54 M55 M56;
     M61 M62 M63 M64 M65 M66];

 % Find the elements of the centripital coupling matrix.
C11 = Tau_l1 - (M(1,:) * ddq.' + G11);
C21 = Tau_l2 - (M(2,:) * ddq.' + G21);
C31 = Tau_l3 - (M(3,:) * ddq.' + G31);
C41 = Tau_l4 - (M(4,:) * ddq.' + G41);
C51 = Tau_l5 - (M(5,:) * ddq.' + G51);
C61 = Tau_l6 - (M(6,:) * ddq.' + G61);
 
G = simplify([G11; G21; G31; G41; G51; G61]);
C = simplify([C11; C21; C31; C41; C51; C61]);

% Print results.
M
C
G

function A = get_A(L, q, dq, ddq, i)
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t)

A = diff(L, dq(i));
At = subs(A, [q dq], ...
    [th1 th2 th3 th4 th5 th6 diff(th1(t), t), diff(th2(t), t), ...
    diff(th3(t), t), diff(th4(t), t), diff(th5(t), t), diff(th6(t), t)]);
dAt = diff(At, t);
A = subs(dAt, [th1 th2 th3 th4 th5 th6 diff(th1(t), t), diff(th2(t), t), ...
    diff(th3(t), t), diff(th4(t), t), diff(th5(t), t), diff(th6(t), t)...
    diff(th1(t),t,t), diff(th2(t),t,t), diff(th3(t),t,t), diff(th4(t),t,t), ...
    diff(th5(t),t,t), diff(th6(t),t,t)], [q dq ddq]);
end

function Kt = trans_kinetic(dth, Jv, m)
Kt = 0.5*dth*Jv.'*m*Jv*dth.';
end

function Kw = rot_kinetic(dth, Jw, R, I)
Kw = 0.5*dth*Jw.'*R*I*R.'*Jw*dth.';
end

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