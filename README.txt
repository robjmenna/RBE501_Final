# RBE501_Final
RBE501 Final Project

Authors: Robert Menna, Michael DeMalia, Gregory Kashmanian

This matlab project lays out the dynamic model for the Fanuc articulated arm robot. The entry point of the application is main.m. 
This script creates six different trajectories based on the initial and final joint positions defined in q0 and qf. From there, 
the script calls generate_dynamic_model.m which will generate a symbolic model of the robot based on the characteristics in the 
data sheet and the payload passed in to the functions inputs. This model can then be passed to generate_torque.m along with the
previously generated pathing, which will generate the required torque at each joint to achieve the desired position. Each joint 
torque is then plotted at the end of the script.

This script normally takes a very long time to run (~20 min) in its default configuration. In order to shorten the run time, the
rotational kinetic energy can be removed from the model and speed up the execution time significantly. Simply, comment out lines 
64,65, and 66 from generate_dynamic_model.m and then remove the rotational terms from the total kinetic energy in line 75 of the 
same file.