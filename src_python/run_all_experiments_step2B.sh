#!/usr/bin/env bash
# Winning model from Step 2A was reduced 3.0, qval 0.05, remove N5
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --nosv --qval 0.05 --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --nocna --qval 0.05 --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --nosv --nocna --qval 0.05 --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --nofocal --qval 0.05 --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --noarms --qval 0.05 --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --nomuts --qval 0.05 --remove_largest_n 5 &