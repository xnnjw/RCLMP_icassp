# Robust Covairance-Learning based Matching Pursuit (RCL-MP)

## Description
This MATLAB code simulates activity detection in massive MIMO systems using different algorithms (CWO, CL-OMP, and RCL-MP) under various scenarios. The script performs Monte Carlo simulations to evaluate the performance metrics (PMD and Per) over different numbers of antennas (M) and pilots (L) for varying numbers of active users (K). The results are visualized in plots.

## Files
- **`functions/`**: Contains helper functions used in the main script.
- **`utils/`**: Utility functions that assist with data processing and simulation setup.
- **`main_script.m`**: The main script that runs the activity detection simulations and generates plots.

## Usage
1. Ensure that the necessary functions and utils are in the respective folders (`functions/` and `utils/`).
2. Run the `main_script.m` to perform the simulations and plot the results.
3. The script will automatically simulate and plot the following:
   - PMD and Per vs. the number of antennas (M) for different values of K.
   - PMD and Per vs. the number of pilots (L) for different values of K.

## General Settings
- **MC_iters**: Number of Monte Carlo iterations (default is 1).
- **P**: Total power scale of the input signal (default is 1).
- **N**: Number of Machine Type Devices (MTDs) (default is 1024).
- **fading**: Fading type for the channel, set to 'uniform' with limits from -15 to 0 dB.
- **pilot**: Pilot sequence type, set to 'bernoulli'.

## Simulation Procedure
### 1. Over M (Antenna Count)
- Number of pilots (L) is fixed at 32.
- Vary the number of antennas (M) over the range [8, 16, 24, 32].
- Evaluate PMD and Per for different numbers of active users (K).

### 2. Over L (Pilot Count)
- Number of antennas (M) is fixed at 32.
- Vary the number of pilots (L) over the range [8, 16, 24, 32].
- Evaluate PMD and Per for different numbers of active users (K).

## Plotting
- The results are plotted with logarithmic scales for PMD.
- Plots display PMD vs. M (for different K) and PMD vs. L (for different K).
- Legends indicate the algorithms used in the comparison.

## Outputs
- **Tensor2Plot_pmd_M**: Performance matrix for PMD over different values of M.
- **Tensor2Plot_pmd_L**: Performance matrix for PMD over different values of L.
- Figures showing the performance metrics.

## Example Command to Run
Simply run the main script: 
```matlab
>> main_script
```
This will execute the full simulation and generate the required plots.
