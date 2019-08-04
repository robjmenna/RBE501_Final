function T0f = forward_kinematics(q)

% Find the forward kinematics.
d = [330 0 0 335 0 80];
theta = [q(1) q(2)+pi/2 q(3) q(4) q(5) q(6)];
a = [50 330 35 0 0 0];
alpha = [pi/2 0 pi/2 -pi/2 pi/2 0];

T01 = DH(alpha(1), a(1), d(1), theta(1));
T12 = DH(alpha(2), a(2), d(2), theta(2));
T23 = DH(alpha(3), a(3), d(3), theta(3));
T34 = DH(alpha(4), a(4), d(4), theta(4));
T45 = DH(alpha(5), a(5), d(5), theta(5));
T5f = DH(alpha(6), a(6), d(6), theta(6));

T02 = T01 * T12;
T03 = T02 * T23;
T04 = T03 * T34;
T05 = T04 * T45;
T0f = T05 * T5f;

end