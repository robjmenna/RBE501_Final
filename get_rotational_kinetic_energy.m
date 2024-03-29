%%
% Find the translational kinetic energy of a link on the robot.

function Kw = get_rotational_kinetic_energy(dth, Jw, R, I)
Kw = 0.5*dth*Jw.'*R*I*R.'*Jw*dth.';
end