#!/usr/bin/env bash
# Winning model from Step 2A was reduced 3.2, remove N5
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nosv &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nocna &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nosv --nocna &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nofocal &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --noarms &
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --nomuts &