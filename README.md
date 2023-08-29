# Redundancy_of_manipulator

[Code type]
Matlab

[Development period]
2023.05 ~ 06

[Purpose]
Do multi secondary tasks such as 'Avoiding joint limit', 'Obstacle avoidance' using redundancy of manipulator
To do that, I use 'Pseudo Inverse Kinematics(Pseudo-IK)', 'Weighted Generalized IK', 'Gradient Projection Method(GPM)'.
And I simulated in Matlab and checked if secondary tasks were alse satisfied in environments including stationary and dynamic obstacles. 

[Result]
GPM was the best among them but it still had limitation that it's quite tricky to set cost function and parameters. 
That means GPM was sensitive for tuning them.
