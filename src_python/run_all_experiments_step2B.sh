#!/usr/bin/env bash
# Winning model from Step 2A was reduced 3.2, remove N5
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nosv --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nocna --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nosv --nocna --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nofocal --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --noarms --remove_largest_n 5 &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nomuts --remove_largest_n 5 &