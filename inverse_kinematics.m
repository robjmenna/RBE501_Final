%%
% Get the joint positions of the robot given a coordinate in the world frame.

function theta_matrix = inverse_kinematics(xc, yc, zc)
%RBE-501
%Term Project
%Mike DeMalia, Gregory Kashmanian and Robert Menna
%Inverse Kinematics

%Defining known variables that form the T60 matrix
%yc=0;
%xc=465;
%zc=695;

r11=1;
r21=0;
r31=0;
r12=0;
r22=-1;
r32=0;
r13=0;
r23=0;
r33=1;

%Defining the T60 matrix
T60=[r11,r21,r31,xc;r12,r22,r32,yc;r13,r23,r33,zc;0,0,0,1];

%Defining the resolved robotic system link lengths
L1 = sqrt((457-127)^2 + 50^2);
L2 = 330;
L3 = sqrt(335^2 + 35^2);

%Location of the end effector from the T60 Matrix
EE_final = T60(1:3,4);

%Rotational 3X3 from T60 matrix
R60 = T60(1:3,1:3);

%End effector offset
EE_Offset = 80;

%Location on the robotic system prior to the end effector
EE_start = EE_final - EE_Offset * R60 * [1 ; 0 ; 0];

%Calculation of theta 1
th1 = atan2(EE_final(2), EE_final(1));

%Positional data at the first joint - Note specific parameters for CR7
%added
P_1 = [50 * cos(th1) ; 50 * sin(th1) ; 330];

%Positional difference from end effector to first joint intersection. 
P_3 = EE_start - P_1;
R_norm = norm([P_3(1), P_3(2), P_3(3)]);

%Calculation of theta 2
P = atan2(P_3(3), norm([P_3(1),P_3(2)]));
G = (acos((R_norm^2 + L2^2 - L3^2) / (2*R_norm*L2)));
th2 = G + P - pi/2;

%Calculation of theta 3 - Note specific parameters for CR7
Alpha = (acos((L2^2 + L3^2 - R_norm^2) / (2*L2*L3)));
Th3_O = atan2(35, 335)
th3 = Alpha - pi/2 - Th3_O;

%Defining parameters needed for theta 4, theta 5 and theta 6
s3=sin(th3);
c3=cos(th3);
c1=cos(th1);
c23=cos(th2+th3);
s1=sin(th1);
s23=sin(th2+th3);

%Calculation of theta 4, 5 and 6
th4=(pi/2)-(atan2(-c1*s23*r13-s1*s23*r23+c23*r33,c1*c23*r13+s1*c23*r23+s23*r33));
th5=(pi/2)-(atan2(sqrt(1-(s1*r13-c1*r23)^2), s1*r13-c1*r23));
c4=cos(th4);
s4=sin(th4);
th6=(pi/2)-(atan2(s1*r12-c1*r22, -s1*r11+c1*r21));

%Overall theta calculation
theta_matrix=(180/pi)*[th1;th2;th3;th4;th5;th6];
end