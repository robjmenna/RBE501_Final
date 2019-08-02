%%
% Get the torque required for each joint to follow a given trajectory.

function Tau = generate_torque(path_1, path_2, path_3, path_4, path_5, path_6, q, dq, ddq, model)
last_index = size(path_1,2);
Tau = zeros(6,last_index);
for i=1:last_index
    fprintf(['Processing torque trajectory. ' num2str(i) '/' num2str(last_index) '\n']);
    
    Tau(1:6,i) = vpa(subs(model,...
        [q dq ddq],...
        [path_1(1,i) path_2(1,i) path_3(1,i) path_4(1,i) path_5(1,i) path_6(1,i)...
        path_1(2,i) path_2(2,i) path_3(2,i) path_4(2,i) path_5(2,i) path_6(2,i)...
        path_1(3,i) path_2(3,i) path_3(3,i) path_4(3,i) path_5(3,i) path_6(3,i)]),3);
end
end