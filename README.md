# Non-linear-predictive-control-NMPC

The main objective is to set : a non linear predictive controller.


Description of the process to be controlled :
The process is a pendulum. The torque U N·m will be manipulated in 
order to control the angle THETA rad. See figure.jpg to understand
the theoritical work given by Prof. André Desbiens.  






Specifications desired : 
- The parameters Hp, Hc and Lambda are adjustable.
- Future set points are assumed to be equal to the current set point.
- The motor cannot supply torques smaller than -1.5 N·m or larger than 1.5 N·m.
- No static error following a step disturbance.
- The observer : an extended Kalman filter.
- Predictions based on the nonlinear model.
