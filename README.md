ML-model-for-Pseudo-Two-Dimensional-P2D-lithium-ion-battery-prediction
1. Workflow
To streamline parameter estimation, I implemented a machine learning (ML) framework that:
1. Generates synthetic data by calling the MATLAB P2D model at different parameter sets.
Trains a deep neural network (surrogate) to learn the mapping.
Identifies parameters that best match a provided experimental dataset, using the
trained surrogate in an optimization loop.
Validates those parameters by plugging them back into the full P2D model and
comparing the simulated voltage with experimental measurements.
This aligns with the Project Objectives:
Develop a Machine Learning Framework for parameter estimation.
Parameter Identification of key P2D parameters that are difficult to measure directly.
Model Validation by comparing predicted voltages to real experimental data.
Optimization to ensure physically plausible parameter values and avoid overfitting.
2. Machine Learning Model Development
2.1 Generating Synthetic Data from the P2D Solver
A central challenge is to generate enough labeled data for training. I directly called the
MATLAB P2D model (MAIN_I_ROM_V3_1_1_PE.m) to obtain the voltage curve for
randomly sampled parameters. Key steps:
2.2 Skipping Invalid (Complex) Solutions
Early tests revealed that some parameter draws drove the PDE solver into unphysical or
numerically unstable regimes, yielding complex voltages. Rather than collecting them and
discarding later, I immediately skip such parameter sets:
This ensures our final dataset is purely real and physically valid.
2.3 Deep Neural Network (Surrogate)
The architecture features multiple dense layers (256–128–64 neurons), batch normalization
for training stability, dropout to reduce overfitting, and EarlyStopping / ReduceLROnPlateau
callbacks:
I then trained the network on our synthetic data, achieving a surrogate that approximates the
PDE solver’s input–output relationship across the parameter/time domain.
3. Parameter Identification & Optimization
3.1 Objective Function
I utilized the trained surrogate in a least-squares sense to match the experimental voltage
data. For each candidate parameter vector:
1. Evaluate the surrogate across all experimental time points.
3.2 Minimization
Using scipy.optimize.minimize with bounds ensures parameters remain within
The result yields best
_params that minimize the SSE to your experimental dataset.
4. Model Validation
4.1 Surrogate vs. Experiment
First, I compare the surrogate predictions at best_params with the experimental data.
This helps confirm that the neural network’s “fast approximation” to the PDE solver matches
the measured data well at least in a local SSE sense.
4.2 True P2D vs. Experiment
To ensure the true PDE model (not just the surrogate) is accurate with our fitted parameters,
I do a final call:
t_sim_final, V_sim_final = simulate_p2d_model(eng, best_params,
data_file, crate_value)
plt.plot(t_sim_final, V_sim_final,
'g-
'
, label=
'Real P2D w/ Best
Params')
Comparing that final curve to the experiment confirms whether the identified parameters are
physically consistent.
Important Note :
I have currently chosen to generate only 50 parameter sets from the P2D
model, primarily because each simulation can be time-consuming and certain
parameter draws occasionally yield complex (invalid) solutions. However,
increasing the number of valid parameter sets would naturally improve the
surrogate’s accuracy and generalizability, as a larger dataset provides better
coverage of the parameter space and more robust training for the neural
network.

<img width="555" alt="Screenshot 2025-04-02 at 9 22 52 PM" src="https://github.com/user-attachments/assets/f1143ee2-6246-4edf-a8a9-9a9c043f7bfe" />
<img width="572" alt="Screenshot 2025-04-02 at 9 23 35 PM" src="https://github.com/user-attachments/assets/1b226b45-fc02-4f96-bd6d-e8301fbc6837" />
<img width="554" alt="Screenshot 2025-04-02 at 9 23 51 PM" src="https://github.com/user-attachments/assets/b9122372-7e8f-4142-99d4-6b8b52f12511" />



