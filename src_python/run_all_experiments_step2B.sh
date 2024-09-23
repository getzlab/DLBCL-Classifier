#!/usr/bin/env bash
# Winning model from Step 2A was reduced 3.2, remove N5
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --nosv &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --nocna &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --nosv --nocna &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --nofocal &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --noarms &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --nomuts &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 --remove_largest_n 5 --noarms --nosv