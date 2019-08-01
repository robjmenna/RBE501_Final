% Calculate an element of the angular component of the Jacobian.
function z = get_Z(T)
    R = T(1:3,1:3);
    z = R * [0;0;1];
end