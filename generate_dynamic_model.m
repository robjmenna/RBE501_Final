function Tau = generate_dynamic_model(q, dq, ddq, payload)

syms a1 a2 a3 d1 d2 d3
fprintf('Starting torque controller\n');

mt = 53;
m1 = mt*0.2869;
m2 = mt*0.3763;
m4 = mt*(0.1713 + 0.1363);
m6 = mt*0.0283;

I1 = 0.47;
I2 = 0.47;
I3 = 0.15;
mp = payload;
g = 9.8;

d = [d1 0 0 d2 0 d3];
theta = [q(1) q(2)+pi/2 q(3) q(4) q(5) q(6)];
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

J_l1 = jacobian(T01(1:3,4), q(1));
J_l2 = jacobian(T02(1:3,4),[q(1) q(2)]);
J_l4t = jacobian(T04(1:3,4),[q(1) q(2) q(3) q(4)]);
J_l6t = jacobian(T0f(1:3,4),[q(1) q(2) q(3) q(4) q(5) q(6)]);
J_l4w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04)];
J_l5w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04) get_Z(T05)];
J_l6w = [get_Z(T01) get_Z(T02) get_Z(T03) get_Z(T04) get_Z(T05) get_Z(T0f)];

v_m1 = J_l1 * dq(1);
v_m2 = J_l2 * [dq(1) ; dq(2)];
v_m4 = J_l4t * [dq(1) ; dq(2) ; dq(3); dq(4)];
v_m6 = J_l6t * [dq(1) ; dq(2) ; dq(3); dq(4); dq(5); dq(6)];

fprintf('Finding K and P\n');
K1 = simplify(subs(0.5 * m1 * (v_m1.' * v_m1), [a1 d1], [a1/2 d1/2]));
K2 = subs(0.5 * m2 * (v_m2.' * v_m2), a2, a2/2);
K4 = subs(0.5 * m4 * (v_m4.' * v_m4), [a3 d2], [a3/2 d2/2]);
K6 = subs(0.5 * m6 * (v_m6.' * v_m6), d3, d3/2);
K_p = 0.5 * mp * (v_m6.' * v_m6);
%K4r = get_rotational_kinetic_energy([dq(1) dq(2) dq(3) dq(4)],J_l4w,T04(1:3,1:3),I4);
%K5r = get_rotational_kinetic_energy([dq(1) dq(2) dq(3) dq(4) dq(5)],J_l5w,T05(1:3,1:3),I5);
%K6r = get_rotational_kinetic_energy([dq(1) dq(2) dq(3) dq(4) dq(5) dq(6)],J_l6w,T0f(1:3,1:3),I6);

P1 = subs(m1 * g * T01(3,4), [a1 d1], [a1/2 d1/2]);
P2 = subs(m2 * g * T02(3,4), a2, a2/2);
P4 = subs(m4 * g * T04(3,4), [a3 d2], [a3/2 d2/2]);
P6 = subs(m6 * g * T0f(3,4), d3, d3/2);
Pp = mp * g * T0f(3,4);

fprintf('Finding L\n');
K = K1 + K2 + K4 + K6 + K_p; %+ K4r + K5r + K6r;
P = P1 + P2 + P4 + P6 + Pp;

L = K - P;
A1 = get_A(L, q, dq, ddq, 1);
A2 = get_A(L, q, dq, ddq, 2);
A3 = get_A(L, q, dq, ddq, 3);
A4 = get_A(L, q, dq, ddq, 4);
A5 = get_A(L, q, dq, ddq, 5);
A6 = get_A(L, q, dq, ddq, 6);

B1 = diff(L,q(1));
B2 = diff(L,q(2));
B3 = diff(L,q(3));
B4 = diff(L,q(4));
B5 = diff(L,q(5));
B6 = diff(L,q(6));

Tau_l1 = A1 - B1;
Tau_l2 = A2 - B2;
Tau_l3 = A3 - B3;
Tau_l4 = A4 - B4;
Tau_l5 = A5 - B5;
Tau_l6 = A6 - B6;

Tau = vpa(subs([Tau_l1;Tau_l2;Tau_l3;Tau_l4;Tau_l5;Tau_l6;], [a1 a2 a3 d1 d2 d3], [50 330 35 330 335 80]),3);