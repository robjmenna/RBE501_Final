%%
% get the first term of the Lagrangian for a given joint in the robot.

function A = get_A(L, q, dq, ddq, joint_num)
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t);

A = diff(L, dq(joint_num));
At = subs(A,[q dq], [th1 th2 th3 th4 th5 th6...
    diff(th1(t),t) diff(th2(t),t) diff(th3(t), t) diff(th4(t),t) diff(th5(t),t) diff(th6(t), t)]);

dAt = diff(At, t);
A = subs(dAt,[th1 th2 th3 th4 th5 th6...
    diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) diff(th4(t),t) diff(th5(t),t) diff(th6(t),t)...
    diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t) diff(th4(t),t,t) diff(th5(t),t,t) diff(th6(t), t,t)],...
    [q dq ddq]);
end