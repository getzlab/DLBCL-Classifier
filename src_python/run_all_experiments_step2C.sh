#!/usr/bin/env bash
# Winning model from Step 2B was same as 2A

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --remove_largest_n 5 --ploidy --savemodels &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --remove_largest_n 5 --coo --savemodels &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.2 --remove_largest_n 5 --coo --ploidy --savemodels &
