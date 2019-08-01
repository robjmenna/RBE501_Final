function Kt = get_translational_kinetic_energy(dth, Jv, m)
Kt = 0.5*dth*Jv.'*m*Jv*dth.';
end