import nn
import generate_test_train_sets as gt
import format_data as fd
import numpy as np
import torch
import torch.nn as torch_nn
import random
import pandas as pd
from operator import itemgetter
import sklearn.metrics
import argparse
import copy


def main(seed, folds, no_sv=False, no_cna=False, ploidy=False, no_silent=False, binary_labels=False, full_features=False,
         reduced_version=None, verbose=False, early_stopping=False, use_pca=None, trim_features=None,
         coo=False, subset=False, no_muts=False, generate_validation_sets=False,
         output_filename=None, training_history=False, qval=0.10, no_arms=False, no_focals=False, bp_cutoff=None, keep_svs=True,
         training_file=None, remove_outliers=False, weightl265p=False, remove_largest_n=None, nosvbcl6=False):
    # Set seeds
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

    print('CUDA:' + str(torch.cuda.is_available()))

    datafile = "../data_tables/gsm/DLBCL_Staudt_Shipp_CL.for_classifier_training.classifier_subset.fix_sv.fix_ploidy.17-Aug-2022.txt"
    targetfile = "../data_tables/confidence_tables/baseline_probabilities.connectivity_based.sensitivity_power2.Aug_17_2022.tsv"

    train_samples = list(pd.read_csv(training_file, sep='\t', header=None, index_col=0).index)

    data, targets = fd.format_inputs(datafile, targetfile, train_samples,
                                     no_sv=no_sv,
                                     no_cna=no_cna,
                                     ploidy=ploidy,
                                     no_silent=no_silent,
                                     reduced_version=reduced_version,
                                     use_pca=use_pca,
                                     coo=coo,
                                     no_muts=no_muts,
                                     qval=qval,
                                     no_arms=no_arms,
                                     no_focals=no_focals,
                                     bp_cutoff=bp_cutoff,
                                     remove_outliers=remove_outliers,
                                     weightl265p=weightl265p,
                                     remove_largest_n=remove_largest_n,
                                     nosvbcl6=nosvbcl6)

    training_sets, validation_sets = gt.generate_kfold_frames(targetfile, folds, list(data.index), seed)

    if generate_validation_sets:
        dirpath = '../all_validation_sets/' + output_filename
        print('Generating validation set for seed:', seed, 'model:', output_filename)

        for i in range(0, len(validation_sets)):
            fn = dirpath+"/"+output_filename+"_"+str(seed)+"_"+str(i)
            with open(fn, 'w+') as f:
                for ele in validation_sets[i]:
                    f.write(ele+'\n')
        print('Sets for seed ' + str(seed) + ' generated.')
        return None, None, None, None

    print("Training samples count: ", len(data))

    nets = []
    t_inputs_arr = []
    t_targets_arr = []
    v_inputs_arr = []
    v_targets_arr = []
    n_features = 0

    # Set up the training and validation frames
    for k in range(0, folds):
        train_data, train_targets, validation_data, validation_targets = \
            gt.generate_train_validation_frames(data, targets, training_sets, validation_sets, k)
        # Security check - Exit immediately if a not-train set sample is found in this training script.
        if  [i for i in list(train_data.index) if i not in train_samples] or \
            [i for i in list(validation_data.index) if i not in train_samples] or \
            [i for i in list(train_targets.index) if i not in train_samples] or \
            [i for i in list(validation_targets.index) if i not in train_samples]:
            print('TEST SET SAMPLE FOUND IN TRAINING SET, EXITING.')
            exit()

        train_data = train_data.astype(float)
        train_targets = train_targets.astype(float)
        validation_data = validation_data.astype(float)
        validation_targets = validation_targets.astype(float)

        net = nn.Net(10, train_data.shape[1], 5)
        nets.append(net)
        t_inputs_arr.append(train_data.values)
        t_targets_arr.append(train_targets.values)
        v_inputs_arr.append(validation_data.values)
        v_targets_arr.append(validation_targets.values)
        if verbose:
            print("Fold: " + str(k) + " Training size: " + str(train_data.shape[0]))
        n_features = train_data.shape[1]
    if verbose:
        if no_sv:
            print("Training with no SVs.")
        if no_cna:
            print("Training with no CNAs")
        if ploidy:
            print("Training with ploidy")
        if no_silent:
            print("Training with no silent mutations")
        if binary_labels:
            print("Training with categorical one-hot target vectors")

    print("Training on " + str(n_features) + ' features')

    if not binary_labels:
        lr = 0.05
    else:
        lr = 0.01

    if early_stopping:
        num_epoch = 100
    else:
        num_epoch = 50

    optimizers = []
    for k in range(0, folds):
        optimizer = torch.optim.SGD(nets[k].parameters(), lr=lr, weight_decay=.001)
        optimizers.append(optimizer)

    if binary_labels:
        criterion = torch_nn.CrossEntropyLoss()
    else:
        criterion = torch_nn.MSELoss()

    train_inputs_arr = [[] for _ in range(folds)]
    train_labels_arr = [[] for _ in range(folds)]
    validation_inputs_arr = [[] for _ in range(folds)]
    validation_labels_arr = [[] for _ in range(folds)]

    # Make each input into a tensor based on probabilistic labels or not
    for k in range(0, folds):
        for i in range(0, len(t_inputs_arr[k])):
            if binary_labels:
                tinput = torch.FloatTensor([t_inputs_arr[k][i]]).requires_grad_(True)
                arr = t_targets_arr[k][i]
                ttarget = torch.LongTensor([np.argmax(arr)])
            else:
                tinput = torch.FloatTensor([t_inputs_arr[k][i]]).requires_grad_(True)
                ttarget = torch.FloatTensor(t_targets_arr[k][i])
            train_inputs_arr[k].append(tinput)
            train_labels_arr[k].append(ttarget)
        for i in range(0, len(v_inputs_arr[k])):
            if binary_labels:
                vinput = torch.FloatTensor([v_inputs_arr[k][i]]).requires_grad_(False)
                arr = v_targets_arr[k][i]
                vtarget = torch.LongTensor([np.argmax(arr)]).requires_grad_(False)
            else:
                vinput = torch.FloatTensor([v_inputs_arr[k][i]]).requires_grad_(False)
                vtarget = torch.FloatTensor(v_targets_arr[k][i]).requires_grad_(False)
            validation_inputs_arr[k].append(vinput)
            validation_labels_arr[k].append(vtarget)

    trained_nets = []
    all_predicted_confidences = []
    all_actual_confidences = []
    predictions = []
    patience = 5
    count = 0
    training_losses = [[] for _ in range(5)]
    validation_losses = [[] for _ in range(5)]

    for k in range(0, folds):
        current_min = float('inf')
        net = nets[k]
        optimizer = optimizers[k]
        curr_train_inputs = train_inputs_arr[k]
        curr_train_labels = train_labels_arr[k]
        initavgloss = float('inf')
        initavgloss_validation = float('inf')
        # rollback model allows us to take the best model before the stopping condition was met (patience-1 steps ago)
        rollback_models = [None] * (patience + 1)

        for epoch in range(0, num_epoch):
            count = (count + 1) % patience
            permutation = list(np.random.permutation(len(curr_train_inputs)))
            curr_train_inputs = [curr_train_inputs[i] for i in permutation]
            curr_train_labels = [curr_train_labels[i] for i in permutation]
            avgloss = 0

            # iterate over all current training inputs, no mini batch
            # zero gradient -> forward inputs -> compute loss -> backprop -> optimizer step
            for i in range(0, len(curr_train_inputs)):
                x = curr_train_inputs[i]
                y = curr_train_labels[i]
                optimizer.zero_grad()
                out = net.forward(x)
                loss = criterion(out, y)
                avgloss += loss.item()
                loss.backward()
                optimizer.step()

            rollback_models.pop()
            rollback_models.insert(0, copy.deepcopy(net))

            # determine if early stopping should occur
            if early_stopping:
                # pass through the validation set
                validation_avg_loss = 0
                curr_validation_inputs = validation_inputs_arr[k]
                curr_validation_labels = validation_labels_arr[k]
                unsorted_outs = []

                for i in range(0, len(curr_validation_inputs)):
                    x = curr_validation_inputs[i]
                    y = curr_validation_labels[i]
                    out = net.forward(x)
                    loss = criterion(out, y)
                    validation_avg_loss += loss.item()
                    output_confidence = np.max(out[0].detach().numpy()).item()
                    actual_confidence = np.max(y.numpy())
                    actual_cluster = np.argmax(y).item()
                    if binary_labels:
                        actual_cluster = y.item()
                    pred_cluster = np.argmax(out[0].detach().numpy()).item()
                    currentlabel = validation_sets[k][i]
                    unsorted_outs.append((output_confidence, actual_confidence, actual_cluster, pred_cluster, currentlabel))

                sorted_outs = sorted(unsorted_outs, key=itemgetter(0))
                actuals = np.array([l[2] for l in sorted_outs])
                preds = np.array([l[3] for l in sorted_outs])

                avgloss = avgloss / len(train_inputs_arr[k])
                validation_avg_loss = validation_avg_loss/len(validation_labels_arr[0])

                if epoch == 0:
                    initavgloss = avgloss
                    initavgloss_validation = validation_avg_loss

                accuracy = sum(actuals == preds)/len(actuals)
                idx = int(len(actuals)*0.3)
                topacc = sum(actuals[idx:] == preds[idx:])/len(actuals[idx:])
                bottomacc = sum(actuals[:idx] == preds[:idx])/len(actuals[:idx])

                if epoch % (num_epoch/10) == 0:
                    print(  'Fold [%d/%d] Epoch [%d/%d] Initial Loss: %.4f, Training Loss: %.4f, '
                            'Initial Validation Loss: %.4f, Validation Loss: %.4f,'
                            'V Acc: %0.4f, Top 70th Acc: %0.4f, Bottom 30th Acc: %0.4f' % (k+1, folds,
                                                                                           epoch, num_epoch,
                                                                                           initavgloss, avgloss,
                                                                                           initavgloss_validation, validation_avg_loss,
                                                                                           accuracy,
                                                                                           topacc, bottomacc))
                training_losses[k].append(avgloss)
                validation_losses[k].append(validation_avg_loss)

                # Early stopping checks - if we haven't seen improvement in patience (5) epochs, end training
                # and remove the last 5 epochs from the training and validation evaluation histories
                if validation_avg_loss > current_min and count == 0:
                    print('Stopping early at Epoch (base 0):', epoch,
                          ' Rolling back to Epoch (base 0):', epoch - patience,
                          ' Top Validation Acc:', topacc)
                    for i in range(patience):
                        training_losses[k].pop()
                        validation_losses[k].pop()
                    break
                elif validation_avg_loss < current_min:
                    count = 0
                    current_min = validation_avg_loss
                else:
                    continue

        # training is finished, the best model is at the end of rollback_models. point net variable towards model 6 epochs ago
        net = rollback_models.pop()

        if verbose:
            print('Evaluating validation set')
        curr_validation_inputs = validation_inputs_arr[k]
        curr_validation_labels = validation_labels_arr[k]
        unsorted_outs = []
        for i in range(0, len(curr_validation_inputs)):
            x = curr_validation_inputs[i]
            y = curr_validation_labels[i]
            out = net.forward(x)
            output_confidence = np.max(out[0].detach().numpy()).item()
            actual_confidence = np.max(y.numpy())
            actual_cluster = np.argmax(y).item()
            if binary_labels:
                actual_cluster = y.item()
            pred_cluster = np.argmax(out[0].detach().numpy()).item()
            currentlabel = validation_sets[k][i]
            if currentlabel == 'DLBCL_MAYO_DLBCL_3653':
                print(seed, k, list(out.detach().numpy()))
            unsorted_outs.append((output_confidence, actual_confidence, actual_cluster, pred_cluster, currentlabel,
                                 list(out.detach().numpy())))

        sorted_outs = sorted(unsorted_outs, key=itemgetter(0))
        fold_act_confidences = np.array([l[1] for l in sorted_outs])
        fold_pred_confidences = np.array([l[0] for l in sorted_outs])
        all_predicted_confidences = np.append(all_predicted_confidences, fold_pred_confidences)
        all_actual_confidences = np.append(all_actual_confidences, fold_act_confidences)
        [predictions.append(x) for x in sorted_outs]
        trained_nets.append(net)
    if training_history:
        dirpath = '../model_training_history/' + output_filename
        traininghistorydf = pd.DataFrame(training_losses)
        validationhistorydf = pd.DataFrame(validation_losses)
        traininghistorydf.index = ['TrainingFold' + str(k + 1) for k in range(folds)]
        validationhistorydf.index = ['ValidationFold' + str(k + 1) for k in range(folds)]
        concatdf = pd.concat([traininghistorydf, validationhistorydf])
        concatdf.index.name = 'epoch'
        fn = dirpath + '/seed' + str(seed) + '_' + output_filename + '.tsv'
        concatdf.to_csv(fn, sep='\t')
    return trained_nets, all_actual_confidences, all_predicted_confidences, predictions
