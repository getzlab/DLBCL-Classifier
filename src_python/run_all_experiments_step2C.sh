#!/usr/bin/env bash
# Winning model from Step 2B was same as 2A

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --ploidy --savemodels &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --coo --savemodels &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --coo --ploidy --savemodels &
