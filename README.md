# Salp-Simulations
Package of MatLab code that simulates and animates motion of salp-like swimmers.  These swimmers are chains of individual agents in a drag-dominated fluid environment that are elastically connected and move by exerting thrust into the fluid.

![](https://github.com/OSU-LRAM/Salp-Simulations-LRAM/blob/main/SalpAnimations/1_medLinks_cleanedUp.gif)

# Stability Studies

Studies joint stability of salp at various body stiffnesses.  Vector field represents shape velocity resulting from combo link thrust and joint elasticity.  Colored lines represent sample evolutions of salp shape from various initial conditions.

![](https://github.com/OSU-LRAM/Salp-Simulations/blob/main/SalpAnimations/Stability%20vs%20Stiffness.gif)

# How To Use

Make sure you have the dependencies installed: [sysplotter](https://github.com/OSU-LRAM/GeometricSystemPlotter) and [the MatLab gif Add-On](https://www.mathworks.com/matlabcentral/fileexchange/63239-gif).  Sysplotter needs to be initialized so that this code has access to Sysplotter's kinematics utilities.

For a basic example of a salp simulation run 'runQuickSimulation.m'

This file will process a basic simulation of a salp with listed options for varying simulation parameters

