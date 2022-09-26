#!/usr/bin/env bash
# 1 (Step 1) + 3 (qv) + 3 (remove n) + 9 (qval + remove n) total models = 16

for qv in "0.05" "0.01" "0.001"
do
    nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --qval $qv &
done

for n in "5" "20" "80"
do
        if [[ "$n" == "5" ]]
        then
            nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --remove_largest_n $n &
        else
            nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --remove_largest_n $n &
        fi
done

for qv in "0.05" "0.01" "0.001"
do

    for n in "5" "20" "80"
    do
        nohup python3 experiment_driver_nn.py --numiter 100 --folds 5 --earlystopping --reduced 3.3 --qval $qv --remove_largest_n $n &
    done
done