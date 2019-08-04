# RBE501_Final
RBE501 Final Project

Authors: Robert Menna, Michael DeMalia, Gregory K

This matlab project lays out the dynamic model for the Fanuc articulated arm robot. The entry point of the application is main.m. 
This script creates six different trajectories based on the initial and final joint positions defined in q0 and qf. From there, 
the script calls generate_dynamic_model.m which will generate a symbolic model of the robot based on the characteristics in the 
data sheet and the payload passed in to the functions inputs. This model can then be passed to generate_torque.m along with the
previously generated pathing, which will generate the required torque at each joint to achieve the desired position. Each joint 
torque is then plotted at the end of the script.
