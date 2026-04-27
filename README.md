
Project: Effects of Noise on a Time-Delayed Model of Human Stick Balancing
Author : Md Shafiul Alom Khan
Course : EBME 419 (Spring 2026)



This folder contains the MATLAB code used to produce every figure in the
report (main.tex). The simulation is a pendulum-cart model with a time delay
of tau = 0.23 s, predictor feedback control, a sensory dead zone Pi, and
optionally Gaussian motor noise on the fingertip velocity. The base script
Fig5_milton_simu.m was provided by Milton et al. (J. R. Soc. Interface,
2016); all other scripts in this folder are extensions written for this
project.


All stochastic results use the global seed rng(42), which is set ONCE at
the top of compute_noiseIC.m and noise_DZ2_run.m. 


FILE DIRECTORY


Script                        Description                                       Figures produced
----------------------------  ------------------------------------------------  -----------------------------------------------
Fig5_milton_simu.m            base sim, reproduces Milton et al. Fig 5(c,d)     Fig 2  (O1)    fig_reproduction.png
sensitivity_data.m            base script for one IC sensitivity run            (called by compute_sensitivity.m)
compute_sensitivity.m         drives sensitivity_data.m, computes lambda        Fig 3  (O2)    fig_traj.png
                                                                                Fig 4  (O2)    fig_phase.png
deadzone_IC.m                 base script for one dead-zone run                 (called by compute_deadzone.m)
compute_deadzone.m            drives deadzone_IC.m, sweeps Pi in [0.5, 2.0]     Fig 5  (O3)    deadzone.png
                                                                                Fig 6  (O3)    deadzone_sweep.png
                                                                                Fig 7  (O3)    deadzone_sweepb.png
noise_IC.m                    base script for one stochastic run                (called by compute_noiseIC.m, noise_DZ2_run.m,
                                                                                 visualize_noise.m, noiseDZ2plot.m)
compute_noiseIC.m             Monte Carlo sweep, Pi = 0.8 deg, N = 1000         (no figures; produces noise_IC.mat)
noise_DZ2_run.m               Monte Carlo sweep, Pi = 2.0 deg, N = 1000         (no figures; produces noise_DZ2.mat)
visualize_noise.m             plotting script for the Pi = 0.8 deg sweep        Fig 8  (O4)    fig_noise_early.png
                                                                                Fig 9  (O4)    fig_noise_bt.png
                                                                                Fig 10 (O4)    fig_noise_traj.png
                                                                                Fig 11 (O4)    fig_noise_psd.png
noiseDZ2plot.m                plotting script for the combined comparison       Fig 12 (O3+O4) fig_combined_bt.png
                                                                                Fig 13 (O3+O4) fig_combined_traj.png


Saved data (.mat)
-----------------
sensitivity_data.mat          output of compute_sensitivity.m
deadzone_IC.mat               output of compute_deadzone.m
noise_IC.mat                  output of compute_noiseIC.m
noise_DZ2.mat                 output of noise_DZ2_run.m
