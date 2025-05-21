# Min-Max NMPC for Robust Nonlinear Tracking Control

This repository implements a robust Min-Max Nonlinear Model Predictive Control (NMPC) approach for tracking control of nonlinear systems subject to bounded disturbances. The implementation compares standard NMPC and Min-Max NMPC approaches using a nonlinear double integrator system.

## Background

Model Predictive Control (MPC) is an advanced control technique that solves an optimization problem at each time step to determine the optimal control input. However, standard MPC implementations may not perform well in the presence of disturbances or uncertainties.

Min-Max NMPC addresses this limitation by explicitly considering worst-case disturbance scenarios in the control design. This robust control approach minimizes the maximum cost across multiple disturbance scenarios, leading to improved tracking performance under uncertain conditions.

## Features

- Implementation of both standard NMPC and Min-Max NMPC for comparison
- Nonlinear double integrator system with customizable parameters
- Disturbance scenario generation for robust control evaluation
- Comprehensive performance evaluation metrics:
  - Position and velocity tracking
  - Control input analysis
  - Tracking error comparison
  - Computational efficiency analysis
  - Phase portrait visualization

## Mathematical Formulation

The Min-Max NMPC formulation solves the optimization problem:

```
min max J(x, u, w)
 u   w
```

where:
- J is the cost function
- x is the state trajectory
- u is the control input sequence
- w represents bounded disturbances

The implemented system follows the dynamics:
```
ẋ₁ = x₂
ẋ₂ = f(x₁, x₂) + u + d
```
where f(x₁, x₂) = -0.5 * x₂ - 0.2 * x₂³ - 0.1 * sin(x₁) and d represents the disturbance.

disturbance scenarios setting:

Our strategy is:
It is divided into extreme scenarios and immediate scenarios. When nd=3, only w is the maximum,  minimum,  zero, three extreme cases. when nd is greater than 3, then a random sample is added to it, which is consistent with the requirements of random generation mentioned in the paper, and it can also be proved through experiments that it has a certain anti-interference performance for nonlinear tasks

## Usage

1. Set parameters in the script header:
   - `run_minmax_nmpc`: Set to true to run Min-Max NMPC
   - `run_standard_nmpc`: Set to true to run standard NMPC
   - `save_results`: Set to true to save simulation results

2. Adjust system parameters as needed:
   - Sampling time and simulation duration
   - NMPC prediction horizon
   - Number of disturbance scenarios
   - System constraints and cost function weights

3. Run the script to execute the simulation and visualize results

## Results

The implementation generates five comparative visualizations:
1. Position and velocity tracking comparison
2. Control input comparison
3. Tracking error and cumulative error comparison
4. Computation time comparison
5. Phase portrait comparison


The analysis also prints performance metrics including mean absolute error, maximum absolute error, and average computation time.
![computation_time_analysis](https://github.com/user-attachments/assets/ec4a4aa0-c580-468e-b451-6c0da7c2f390)
![error_distribution](https://github.com/user-attachments/assets/19bd4ab9-7482-4333-82a1-660e432e3557)
![control_analysis](https://github.com/user-attachments/assets/e1110915-d3af-4dde-bc4a-d34a31d017de)
![tracking_comparison](https://github.com/user-attachments/assets/ff236d46-d118-4ee0-8ad1-99dace3d9e26)
![time_varying_analysis](https://github.com/user-attachments/assets/a1ba825e-9320-448a-af7e-0f515d051153)
![phase_portrait_analysis](https://github.com/user-attachments/assets/72ffd720-2b3b-428c-82bb-a5ba8f118c04)
![image](https://github.com/user-attachments/assets/655b099f-dc67-48a1-b206-65c878948bdd)



## Requirements


- MATLAB (developed and tested with MATLAB R2020b or later)
- Optimization Toolbox (for `fmincon`)



## License

[MIT License](LICENSE) 
