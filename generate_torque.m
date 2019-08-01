function Tau = generate_torque(path_1, path_2, path_3, path_4, path_5, path_6, payload)
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6
syms m1 m2 m4 m6 mp I4 I5 I6
syms a1 a2 a3 d1 d2 d3 g

q = [q1 q2 q3 q4 q5 q6];
dq = [dq1 dq2 dq3 dq4 dq5 dq6];
ddq = [ddq1 ddq2 ddq3 ddq4 ddq5 ddq6];

mt = 53;

fprintf('Starting torque controller\n');

d = [d1 0 0 d2 0 d3];
theta = [q1 q2+pi/2 q3 q4 q5 q6];
a = [a1 a2 a3 0 0 0];
alpha = [pi/2 0 pi/2 -pi/2 pi/2 0];

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

J_l1 = jacobian(T01(1:3,4), q1);
J_l2 = jacobian(T02(1:3,4),[q1 q2]);
J_l4t = jacobian(T04(1:3,4),[q1 q2 q3 q4]);
J_l6t = jacobian(T0f(1:3,4),[q1 q2 q3 q4 q5 q6]);
J_l4w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04)];
J_l5w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04) get_Z(T05)];
J_l6w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04) get_Z(T05) get_Z(T0f)];

v_m1 = J_l1 * dq1;
v_m2 = J_l2 * [dq1 ; dq2];
v_m4 = J_l4t * [dq1 ; dq2 ; dq3; dq4];
v_m6 = J_l6t * [dq1 ; dq2 ; dq3; dq4; dq5; dq6];

fprintf('Finding K and P\n');
K1 = simplify(subs(0.5 * m1 * (v_m1.' * v_m1), [a1 d1], [a1/2 d1/2]));
K2 = subs(0.5 * m2 * (v_m2.' * v_m2), a2, a2/2);
K4 = subs(0.5 * m4 * (v_m4.' * v_m4), [a3 d2], [a3/2 d2/2]);
K6 = subs(0.5 * m6 * (v_m6.' * v_m6), d3, d3/2);
K_p = 0.5 * mp * (v_m6.' * v_m6);
K4r = get_rotational_kinetic_energy([dq1 dq2 dq3 dq4],J_l4w,T04(1:3,1:3),I4);
K5r = get_rotational_kinetic_energy([dq1 dq2 dq3 dq4 dq5],J_l5w,T05(1:3,1:3),I5);
K6r = get_rotational_kinetic_energy([dq1 dq2 dq3 dq4 dq5 dq6],J_l6w,T0f(1:3,1:3),I6);

P1 = subs(m1 * g * T01(3,4), [a1 d1], [a1/2 d1/2]);
P2 = subs(m2 * g * T02(3,4), a2, a2/2);
P4 = subs(m4 * g * T04(3,4), [a3 d2], [a3/2 d2/2]);
P6 = subs(m6 * g * T0f(3,4), d3, d3/2);
Pp = mp * g * T0f(3,4);

fprintf('Finding L\n');
K = K1 + K2 + K4 + K6 + K_p + K4r + K5r + K6r;
P = P1 + P2 + P4 + P6 + Pp;

L = K - P;
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

last_index = size(path_1, 2);
Tau_s = [Tau_l1;Tau_l2;Tau_l3;Tau_l4;Tau_l5;Tau_l6;];
Tau = zeros(6,last_index);

for i=1:last_index
    fprintf(['Processing torque trajectory. ' num2str(i) '/' num2str(last_index) '\n']);
    
    Tau(1:6,i) = vpa(subs(Tau_s,...
        [q dq ddq g a1 a2 a3 d1 d2 d3 m1 m2 m4 m6 mp I4 I5 I6],...
        [path_1(1,i) path_2(1,i) path_3(1,i) path_4(1,i) path_5(1,i) path_6(1,i)...
        path_1(2,i) path_2(2,i) path_3(2,i) path_4(2,i) path_5(2,i) path_6(2,i)...
        path_1(3,i) path_2(3,i) path_3(3,i) path_4(3,i) path_5(3,i) path_6(3,i)...
        9.8 50 330 35 330 335 80 (mt*0.2869) (mt*0.3763) (mt*(0.1713 + 0.1363)) (mt*0.0283) payload...
        0.47 0.47 0.15]));
end
end