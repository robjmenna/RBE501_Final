function Tau = generate_torque(path_1, path_2, path_3, path_4, path_5, path_6, Kv, Kp)
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6
syms m1 m2 m4 m6 mp pi_
syms a1 a2 a3 d1 d2 d3 g

mt = 53;
% m1 = mt * 0.2869;
% m2 = mt * 0.3763;
% m4 = mt * (0.1713 + 0.1363);
% m6 = mt * (0.0283);
% mp = 0;

fprintf('Starting torque controller\n');

d = [d1 0 0 d2 0 d3];
theta = [q1 q2+pi_/2 q3 q4 q5 q6];
a = [a1 a2 a3 0 0 0];
alpha = [pi_/2 0 pi_/2 -pi_/2 pi_/2 0];

T01 = DH(alpha(1), a(1), d(1), theta(1));
T12 = DH(alpha(2), a(2), d(2), theta(2));
T23 = DH(alpha(3), a(3), d(3), theta(3));
T34 = DH(alpha(4), a(4), d(4), theta(4));
T45 = DH(alpha(5), a(5), d(5), theta(5));
T5f = DH(alpha(6), a(6), d(6), theta(6));

T02 = T01 * T12;
T03 = T02 * T23;
T04 = T03*T34;
T05 = T04*T45;
T0f = T05*T5f;

J_l1 = jacobian(T01(1:3,4),[q1]);
J_l2 = jacobian(T02(1:3,4),[q1 q2]);
%J_l3 = jacobian(T03(1:3,4),[q1 q2 q3]);
J_l4 = jacobian(T04(1:3,4),[q1 q2 q3 q4]);
%J_l5 = jacobian(T05(1:3,4),[q1 q2 q3 q4 q5]);
J_l6 = jacobian(T0f(1:3,4),[q1 q2 q3 q4 q5 q6]);

v_m1 = J_l1 * dq1;
v_m2 = J_l2 * [dq1 ; dq2];
%v_m3 = J_l3 * [dq1 ; dq2 ; dq3];
v_m4 = J_l4 * [dq1 ; dq2 ; dq3; dq4];
%v_m5 = J_l5 * [dq1 ; dq2 ; dq3; dq4; dq5];
v_m6 = J_l6 * [dq1 ; dq2 ; dq3; dq4; dq5; dq6];

fprintf('Finding K and P\n');
K1 = subs(0.5 * m1 * (v_m1.' * v_m1), [a1 d1], [a1/2 d1/2]);
K2 = subs(0.5 * m2 * (v_m2.' * v_m2), a2, a2/2);
K4 = subs(0.5 * m4 * (v_m4.' * v_m4), [a3 d2], [a3/2 d2/2]);
K6 = subs(0.5 * m6 * (v_m6.' * v_m6), d3, d3/2);
K_p = 0.5 * mp * (v_m6.' * v_m6);

P1 = subs(m1 * g * T01(3,4), [a1 d1], [a1/2 d1/2]);
P2 = subs(m2 * g * T02(3,4), a2, a2/2);
P4 = subs(m4 * g * T04(3,4), [a3 d2], [a3/2 d2/2]);
P6 = subs(m6 * g * T0f(3,4), d3, d3/2);
Pp = mp * g * T0f(3,4);
fprintf('Finding L\n');
K = K1 + K2 + K4 + K6 + K_p;
P = P1 + P2 + P4 + P6 + Pp;

L = K - P;

q = [q1 q2 q3 q4 q5 q6];
dq = [dq1 dq2 dq3 dq4 dq5 dq6];
ddq = [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6];
A1 = get_A(L, q, dq, ddq, 1);
A2 = get_A(L, q, dq, ddq, 2);
A3 = get_A(L, q, dq, ddq, 3);
A4 = get_A(L, q, dq, ddq, 4);
A5 = get_A(L, q, dq, ddq, 5);
A6 = get_A(L, q, dq, ddq, 6);

B1 = diff(L,q1);
B2 = diff(L,q2);
B3 = diff(L,q3);
B4 = diff(L,q4);
B5 = diff(L,q5);
B6 = diff(L,q6);

Tau_l1 = A1 - B1;
Tau_l2 = A2 - B2;
Tau_l3 = A3 - B3;
Tau_l4 = A4 - B4;
Tau_l5 = A5 - B5;
Tau_l6 = A6 - B6;

fprintf('Simplifying Torque equation.\n');
M11 = (Tau_l1 - subs(Tau_l1,ddq1,0)) /ddq1;
M12 = (Tau_l1 - subs(Tau_l1,ddq2,0)) /ddq2;
M13 = (Tau_l1 - subs(Tau_l1,ddq3,0)) /ddq3;
M14 = (Tau_l1 - subs(Tau_l1,ddq4,0)) /ddq4;
M15 = (Tau_l1 - subs(Tau_l1,ddq5,0)) /ddq5;
M16 = (Tau_l1 - subs(Tau_l1,ddq6,0)) /ddq6;
fprintf('Done M row 1.\n');

M21 = (Tau_l2 - subs(Tau_l2,ddq1,0)) /ddq1;
M22 = (Tau_l2 - subs(Tau_l2,ddq2,0)) /ddq2;
M23 = (Tau_l2 - subs(Tau_l2,ddq3,0)) /ddq3;
M24 = (Tau_l2 - subs(Tau_l2,ddq4,0)) /ddq4;
M25 = (Tau_l2 - subs(Tau_l2,ddq5,0)) /ddq5;
M26 = (Tau_l2 - subs(Tau_l2,ddq6,0)) /ddq6;
fprintf('Done M row 2\n.');

M31 = (Tau_l3 - subs(Tau_l3,ddq1,0)) /ddq1;
M32 = (Tau_l3 - subs(Tau_l3,ddq2,0)) /ddq2;
M33 = (Tau_l3 - subs(Tau_l3,ddq3,0)) /ddq3;
M34 = (Tau_l3 - subs(Tau_l3,ddq4,0)) /ddq4;
M35 = (Tau_l3 - subs(Tau_l3,ddq5,0)) /ddq5;
M36 = (Tau_l3 - subs(Tau_l3,ddq6,0)) /ddq6;
fprintf('Done M row 3.\n');

M41 = (Tau_l4 - subs(Tau_l4,ddq1,0)) /ddq1;
M42 = (Tau_l4 - subs(Tau_l4,ddq2,0)) /ddq2;
M43 = (Tau_l4 - subs(Tau_l4,ddq3,0)) /ddq3;
M44 = (Tau_l4 - subs(Tau_l4,ddq4,0)) /ddq4;
M45 = (Tau_l4 - subs(Tau_l4,ddq5,0)) /ddq5;
M46 = (Tau_l4 - subs(Tau_l4,ddq6,0)) /ddq6;
fprintf('Done M row 4.\n');

M51 = (Tau_l5 - subs(Tau_l5,ddq1,0)) /ddq1;
M52 = (Tau_l5 - subs(Tau_l5,ddq2,0)) /ddq2;
M53 = (Tau_l5 - subs(Tau_l5,ddq3,0)) /ddq3;
M54 = (Tau_l5 - subs(Tau_l5,ddq4,0)) /ddq4;
M55 = (Tau_l5 - subs(Tau_l5,ddq5,0)) /ddq5;
M56 = (Tau_l5 - subs(Tau_l5,ddq6,0)) /ddq6;
fprintf('Done M row 5.\n');

M61 = (Tau_l6 - subs(Tau_l6,ddq1,0)) /ddq1;
M62 = (Tau_l6 - subs(Tau_l6,ddq2,0)) /ddq2;
M63 = (Tau_l6 - subs(Tau_l6,ddq3,0)) /ddq3;
M64 = (Tau_l6 - subs(Tau_l6,ddq4,0)) /ddq4;
M65 = (Tau_l6 - subs(Tau_l6,ddq5,0)) /ddq5;
M66 = (Tau_l6 - subs(Tau_l6,ddq6,0)) /ddq6;
fprintf('Done M row 6.\n');

G11 = subs(Tau_l1, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
G21 = subs(Tau_l2, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
G31 = subs(Tau_l3, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
G41 = subs(Tau_l4, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
G51 = subs(Tau_l5, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
G61 = subs(Tau_l6, [dq1 dq2 dq3 dq4 dq5 dq6 ddq1 ddq2 ddq3 ddq4 ddq5 ddq6], [0 0 0 0 0 0 0 0 0 0 0 0]);
fprintf('Done G.\n');

M = [M11 M12 M13 M14 M15 M16;
M21 M22 M23 M24 M25 M26;
M31 M32 M33 M34 M35 M36;
M41 M42 M43 M44 M45 M46;
M51 M52 M53 M54 M55 M56;
M61 M62 M63 M64 M65 M66];

G = [G11;
G21;
G31;
G41;
G51;
G61];

C11 = Tau_l1 - (M(1,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G11);
C21 = Tau_l2 - (M(2,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G21);
C31 = Tau_l3 - (M(3,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G31);
C41 = Tau_l4 - (M(4,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G41);
C51 = Tau_l5 - (M(5,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G51);
C61 = Tau_l6 - (M(6,:) * [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6].' + G61);
fprintf('Done C.\n');

C = [C11;
C21;
C31;
C41;
C51;
C61];

fprintf('Calculating controller.\n');
last_index = size(path_1, 2);
uno = [path_1(3,last_index); path_2(3,last_index); path_3(3,last_index);...
    path_4(3,last_index); path_5(3,last_index); path_6(3,last_index)] - transpose(ddq);
dos = [path_1(2,last_index); path_2(2,last_index); path_3(2,last_index);...
    path_4(2,last_index); path_5(2,last_index); path_6(2,last_index)] - transpose(dq);
tres = [path_1(1,last_index); path_2(1,last_index); path_3(1,last_index);...
    path_4(1,last_index); path_5(1,last_index); path_6(1,last_index)]- transpose(q);

%a = uno + Kv*dos + Kp*tres;
Tau = zeros(6,last_index);
%quatro = M*a;
controller = M+C+G;
for i=1:last_index
    fprintf(['Processing torque trajectory.' num2str(i) '/' num2str(last_index) '\n']);

    Tau(1:6,i) = subs(controller,...
        [q dq ddq g a1 a2 a3 d1 d2 d3 m1 m2 m4 m6 mp pi_],...
        [path_1(1,i) path_2(1,i) path_3(1,i) path_4(1,i) path_5(1,i) path_6(1,i)...
        path_1(2,i) path_2(2,i) path_3(2,i) path_4(2,i) path_5(2,i) path_6(2,i)...
        path_1(3,i) path_2(3,i) path_3(3,i) path_4(3,i) path_5(3,i) path_6(3,i)...
        9.8 50 330 35 330 335 80 (mt*0.2869) (mt*0.3763) (mt*(0.1713 + 0.1363)) (mt*0.0283) 0 3.14]);
end

end