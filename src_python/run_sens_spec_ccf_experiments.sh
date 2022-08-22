#!/usr/bin/env bash
nohup python3 random_add_in_experiment.py &
nohup python3 random_dropout_experiment.py &
nohup python3 ccf_threshold_experiment.py &