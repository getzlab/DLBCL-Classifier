#!/usr/bin/env bash
nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --fullfeatures &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --fullfeatures &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --fullfeatures &

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.4 &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --reduced 3.4 &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --reduced 3.4 &

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --pca 40 &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --pca 40 &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --pca 40 &

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --pca 30 &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --pca 30 &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --pca 30 &

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --pca 20 &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --pca 20 &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --pca 20 &

nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --pca 2 &
nohup python3 experiment_driver_rf.py --numiter 100 --folds 5 --pca 2 &
nohup python3 experiment_driver_mnb.py --numiter 100 --folds 5 --pca 2 &
