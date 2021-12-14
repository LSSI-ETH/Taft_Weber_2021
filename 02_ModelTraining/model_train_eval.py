#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sept 18 2021

@author: bgao

Code for Taft J
Train RF and RNN models on binding/non-binding RBD sequences generated through
yeast display selections. 


"""

import os
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, matthews_corrcoef, f1_score
from sklearn.utils import shuffle
from tensorflow import keras
from tensorflow.keras import layers
import h5py

import seaborn as sns
import matplotlib.pyplot as plt

import joblib

os.environ['KMP_DUPLICATE_LIB_OK']='True'

random_state = random.randint(1,420)
print("Random State for this run is: "+str(random_state))

from helpers import encode_onehot_padded, flatten_matrix, levenshtein_RBD, anti_merge, make_dist_df, clean_df, concat_dfs, balance_edit_distance, one_hot, RNN_encode_data, RNN_model_test_on_distances, RF_model_test_on_distances


"""
Step 1: Create RNN Functions
"""
def fit_RNN(X_train, y_train, RNN_name):
    # define the model
    model = keras.Sequential()
    model.add(layers.LSTM(units = 80, return_sequences = True, dropout = 0.2, input_shape = (X_train.shape[1],X_train.shape[2])))
    model.add(layers.LSTM(units = 80, return_sequences = True, dropout = 0.2))
    model.add(layers.LSTM(units = 80, return_sequences = False, dropout = 0.2))
    model.add(layers.Dense(units = 50, activation = 'relu'))
    model.add(layers.Dense(units = 1, activation = 'sigmoid'))
    
    model.compile(optimizer = keras.optimizers.Adam(), loss = 'binary_crossentropy',
                  metrics = ['accuracy', keras.metrics.AUC(curve = 'PR')])

    history_full = model.fit(X_train, y_train, batch_size = 32, epochs = 50, verbose = 2, validation_split = 0.1)
    
    plt.plot(history_full.history['accuracy'])
    plt.plot(history_full.history['val_accuracy'])
    plt.title(RNN_name+'_RNN_All_3LSTM_model accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.show()
    plt.savefig(RNN_name+'_RNN_All_3LSTM_Accuracy')
    plt.close()
    
    # summarize history for loss
    plt.plot(history_full.history['loss'])
    plt.plot(history_full.history['val_loss'])
    plt.title(RNN_name+'_RNN_All_3LSTM_model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.show()
    plt.savefig(RNN_name+'_RNN_All_3LSTM_Loss')
    plt.close()
    
    # Save the model
    model.save(RNN_name+"_RNN_All_3LSTM.h5")
    print(RNN_name+' RNN Model has been trained!')
    return model 

def train_RNN(mAb_df, mAb_name, random_state, same_state = True):
    # Create RNN name
    RNN_name = mAb_name+'_RNN'
    # Set random state if we want to shuffle the dataset splits
    if same_state == True:
        rand_state = random_state
    elif same_state == False:
        rand_state = random.randint(1,420)
    
    print("Random State for "+RNN_name+" is: "+str(rand_state))

    # Balance df and create train/test set
    distances = list(range(1,15))
    mAb_df_final = balance_edit_distance(mAb_df, distances, RNN_name, rand_state)
    
    # split into train/test sets
    train_df, test_df = train_test_split(mAb_df_final, test_size = 0.1, random_state = rand_state)
    train_df.to_csv(RNN_name+"_train_data.csv")
    test_df.to_csv(RNN_name+"_test_data.csv")
    
    # Create label from datasets
    y_train_RNN = train_df.Label
    y_test_RNN = test_df.Label
        
    # encode the training data
    X_train_RNN = encode_onehot_padded(train_df.junction_aa)
    X_test_RNN = encode_onehot_padded(test_df.junction_aa)
    
    # Train the model
    trained_RNN = fit_RNN(X_train_RNN, y_train_RNN, RNN_name)

    # Evaluate the model
    y_pred = trained_RNN.predict(X_test_RNN)
    y_pred[y_pred <= 0.5] = 0
    y_pred[y_pred > 0.5] = 1
    model_acc = accuracy_score(y_test_RNN, y_pred)
    model_f1 = f1_score(y_test_RNN, y_pred)
    
    # predict on accuracies
    distances = list(range(1,11))
    model_performance = RNN_model_test_on_distances(distances, test_df, trained_RNN)
    model_performance.to_csv(RNN_name+"_RNN_All_3LSTM_Scores.csv")

    # Create graph of accuracies
    graph_labels = list(map(str, distances))
    graph_positive = model_performance['Pos_Accuracies'].squeeze()
    graph_negative = model_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,ax = plt.subplots()
    rects1 = ax.bar(a-width/2, graph_positive, width, label = 'Positive', color = 'orange')
    rects2 = ax.bar(a+width/2, graph_negative, width, label = 'Negative', color = 'blue')

    ax.set_ylabel('Accuracies')
    ax.set_title(RNN_name + ' RNN Accuracy on Test Sequences by Edit Distance')
    ax.set_xticks(a)
    ax.set_xticklabels(graph_labels)
    ax.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(RNN_name + '_RNN_Accuracy_by_Distance')
    plt.close()

    # Then predict accuracies on unseen data
    unseen_data = pd.read_table(RNN_name+'_unseen.tsv')
    unseen_performance = RNN_model_test_on_distances(distances, unseen_data, trained_RNN)
    unseen_performance.to_csv(RNN_name+"_RNN_All_3LSTM_Unseen_Scores.csv")


    graph_labels = list(map(str, distances))
    graph_positive2 = unseen_performance['Pos_Accuracies'].squeeze()
    graph_negative2 = unseen_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,dx = plt.subplots()
    rects1 = dx.bar(a-width/2, graph_positive2, width, label = 'Positive', color = 'orange')
    rects2 = dx.bar(a+width/2, graph_negative2, width, label = 'Negative', color = 'blue')

    dx.set_ylabel('Accuracies')
    dx.set_title(RNN_name + ' RNN Unseen Accuracy on Test Sequences by Edit Distance')
    dx.set_xticks(a)
    dx.set_xticklabels(graph_labels)
    dx.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(RNN_name + '_RNN_Unseen_Accuracy_by_Distance')
    plt.close()


    # Finally, make the PR-ROC curve
    from sklearn.metrics import roc_curve, auc
    from sklearn.utils import shuffle
    
    y_test_shuff = shuffle(y_test_RNN)
    
    y_pred_keras = trained_RNN.predict_proba(X_test_RNN).ravel()
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test_RNN, y_pred_keras)
    auc_keras = auc(fpr_keras, tpr_keras)

    fpr_keras_shuff, tpr_keras_shuff, thresholds_keras_shuff = roc_curve(y_test_shuff, y_pred_keras)
    auc_keras_shuff = auc(fpr_keras_shuff, tpr_keras_shuff)

    
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_keras, tpr_keras, label=(mAb_name+'_RNN (Full) (area = {:.3f})').format(auc_keras))
    plt.plot(fpr_keras_shuff, tpr_keras_shuff, label=(mAb_name+'_RNN (Shuffled) (area = {:.3f})').format(auc_keras_shuff))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve for '+mAb_name+' Full RNN')
    plt.legend(loc='best')
    plt.show()
    plt.savefig(RNN_name+"_RNN_ROCCurve")
    plt.close()

    # Save model as a pickle    
    print(RNN_name + " RNN Model has been evaluated!")
    return trained_RNN, test_df, model_acc, model_f1

"""
Step 2: Create RF Functions
"""

def RF_Random_Search(X_train, y_train, model, parameters, n_iter, mAb_name):
    from sklearn.model_selection import RandomizedSearchCV
    # Create the classifier
    clf = model
    # Perform Random Search CV
    clf_search = RandomizedSearchCV(clf, param_distributions = parameters,
                                    n_iter = n_iter, cv = 5, verbose = 2,
                                    scoring = 'precision',
                                    return_train_score = True,
                                    n_jobs = -1)
    
    clf_search.fit(X_train, y_train)
    best_parameters = clf_search.best_params_
    best_precision = clf_search.best_score_
    results = clf_search.cv_results_
    best_clf = clf_search.best_estimator_
    
    joblib.dump(best_clf, mAb_name + '_RF_model.gz')
    
    with open(mAb_name + "_GridSearchCV_Best.txt","w") as outfile:
        print(mAb_name + "_Best Parameters: "+str(best_parameters)+'\n',file = outfile)
        print(mAb_name + "_Best Precision: "+str(best_precision)+'\n',file = outfile)
        print(mAb_name + "_Results: "+str(results), file=outfile)
    df = pd.DataFrame(data = results)
    print(df)
    with open(mAb_name + "_GridSearchCV.csv","w") as outfile:
        df.to_csv(outfile)
    return best_clf

def train_RF(mAb_df, parameters, n_iter, mAb_name, random_state, same_state = True):
    # Create RNN name
    RF_name = mAb_name+'_RF'
    # Set random state
    if same_state == True:
        rand_state = random_state
    elif same_state == False:
        rand_state = random.randint(1,420)
        
    print("Random State for "+RF_name+" is: "+str(rand_state))

    # Balance df and create train/test set
    distances = list(range(1,15))
    mAb_df_final = balance_edit_distance(mAb_df, distances, RF_name, rand_state)
    
    # Now split into train/test sets and save
    train_df, test_df = train_test_split(mAb_df_final, test_size = 0.1, random_state = rand_state)
    train_df.to_csv(RF_name+"_train_data.csv")
    test_df.to_csv(RF_name+"_test_data.csv")
    
    # create the model
    model = RandomForestClassifier()
    
    # Create X and y from the split sets
    X_train = train_df.junction_aa
    y_train = train_df.Label
    X_test = test_df.junction_aa
    y_test = test_df.Label
    
    # Create the test_df that we'll use later for distance predictions
    X_test_df = pd.DataFrame()
    X_test_df['junction_aa'] = X_test
    
    # Now encode the training data
    X_train_oh = flatten_matrix(encode_onehot_padded(X_train))
    X_test_oh = flatten_matrix(encode_onehot_padded(X_test))
    
    trained_model = RF_Random_Search(X_train_oh, y_train, model, parameters, n_iter, RF_name)

    # Train the model
    cv_score = cross_val_score(trained_model, X_train_oh, y_train, cv = 5)
    y_pred = trained_model.predict(X_test_oh)
    model_acc = accuracy_score(y_test, y_pred)
    model_f1 = f1_score(y_test, y_pred)
    
    # predict by accuracies
    distances = list(range(1,11))
    test_data = pd.merge(X_test_df, mAb_df, how = "inner", on = 'junction_aa')
    model_performance = RF_model_test_on_distances(distances, test_df, trained_model)

    # Graph accuracies
    graph_labels = list(map(str, distances))
    graph_positive = model_performance['Pos_Accuracies'].squeeze()
    graph_negative = model_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,ax = plt.subplots()
    rects1 = ax.bar(a-width/2, graph_positive, width, label = 'Positive', color = 'orange')
    rects2 = ax.bar(a+width/2, graph_negative, width, label = 'Negative', color = 'blue')

    ax.set_ylabel('Accuracies')
    ax.set_title(RF_name + ' Model Accuracy on Test Sequences by Edit Distance')
    ax.set_xticks(a)
    ax.set_xticklabels(graph_labels)
    ax.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(RF_name + '_Model_Accuracy_by_Distance')
    plt.close()

    model_performance.to_csv(RF_name+'_Model_Accuracies_by_Distance.csv')

    # Then predict accuracies on unseen data
    unseen_data = pd.read_table(RF_name+'_unseen.tsv')
    unseen_performance = RF_model_test_on_distances(distances, unseen_data, trained_model)

    graph_labels = list(map(str, distances))
    graph_positive2 = unseen_performance['Pos_Accuracies'].squeeze()
    graph_negative2 = unseen_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,dx = plt.subplots()
    rects1 = dx.bar(a-width/2, graph_positive2, width, label = 'Positive', color = 'orange')
    rects2 = dx.bar(a+width/2, graph_negative2, width, label = 'Negative', color = 'blue')

    dx.set_ylabel('Accuracies')
    dx.set_title(RF_name + ' Model Unseen Accuracy on Test Sequences by Edit Distance')
    dx.set_xticks(a)
    dx.set_xticklabels(graph_labels)
    dx.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(RF_name + '_Model_Unseen_Accuracy_by_Distance')
    plt.close()
    unseen_performance.to_csv(RF_name+'_Model_Unseen_Accuracies_by_Distance.csv')
    
    # Finally, make the PR-ROC curve
    from sklearn.metrics import roc_curve, auc
    from sklearn.utils import shuffle
    
    y_test_shuff = shuffle(y_test)
    
    y_pred_RF = trained_model.predict_proba(X_test_oh)[:,1].ravel()
    fpr_RF, tpr_RF, thresholds_RF = roc_curve(y_test, y_pred_RF)
    auc_RF = auc(fpr_RF, tpr_RF)
    
    fpr_RF_shuff, tpr_RF_shuff, thresholds_RF_shuff = roc_curve(y_test_shuff, y_pred_RF)
    auc_RF_shuff = auc(fpr_RF_shuff, tpr_RF_shuff)
    
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_RF, tpr_RF, label=(mAb_name+'_RF (Full) (area = {:.3f})').format(auc_RF))
    plt.plot(fpr_RF_shuff, tpr_RF_shuff, label=(mAb_name+'_RF (Shuffled) (area = {:.3f})').format(auc_RF_shuff))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve for Full Models')
    plt.legend(loc='best')
    plt.show()
    plt.savefig(mAb_name+"_RF_ROCCurve")
    plt.close()

    # Save model as a pickle    
    print(mAb_name + " RF Model has been evaluated!")
    return trained_model, test_data, cv_score, model_acc, model_f1


"""
Write a function to train both RNN and RF on the same dataset
"""

def train_both_models(mAb_names, random_state):
    # Define other variables
    RBD_seq = "KNEGFNCYFPLQSYGFQPTNGVGY"
    syn_seq_list = pd.read_csv('46_selected_df_single_mAb.csv')
    
    # Make new function to exclude the 46 sequences from all datasets
    def anti_merge_seqs(df_a, df_b):
        merge_df = pd.merge(df_a, df_b, how = 'outer', on = ['junction_aa'], indicator = 'Ind')
        merge_df_left = merge_df[merge_df['Ind'] == 'left_only']
        return (merge_df_left)

    # Initiate Score Dataframes
    RF_cv_scores = pd.DataFrame()
    RF_accuracy_scores = pd.DataFrame()
    RF_f1_scores = pd.DataFrame()
    
    RNN_accuracy_scores = pd.DataFrame()
    RNN_f1_scores = pd.DataFrame()
    
    # Parameters for RF
    n_estimators = [int(x) for x in np.linspace(start = 50, stop = 500, num  = 10)]
    max_features = ['auto', 'sqrt']
    max_depth = [20,40,60,80,100,150,None]
    min_samples_split = [2, 5, 10]
    min_samples_leaf = [1, 2, 5]

    parameters = {'n_estimators':n_estimators,
                  'max_features':max_features,
                  'max_depth':max_depth,
                  'min_samples_split':min_samples_split,
                  'min_samples_leaf':min_samples_leaf}

    
    for i in range(len(mAb_names)):
        mAb_name = mAb_names[i]
        # Import data for mAb_df
        C2_pos = pd.read_table('2C_'+mAb_name+'_B_df.tsv')
        C2_neg = pd.read_table('2C_'+mAb_name+'_E_df.tsv')
        T2_pos = pd.read_table('2T_'+mAb_name+'_B_df.tsv')
        T2_neg = pd.read_table('2T_'+mAb_name+'_E_df.tsv')

        # Clean the df's
        C2_pos = clean_df(C2_pos, 1, mAb_name, RBD_seq)
        C2_neg = clean_df(C2_neg, 0, mAb_name, RBD_seq)
        T2_pos = clean_df(T2_pos, 1, mAb_name, RBD_seq)
        T2_neg = clean_df(T2_neg, 0, mAb_name, RBD_seq)
        
        # Combine the dfs
        mAb_df = concat_dfs([C2_pos, C2_neg, T2_pos, T2_neg], mAb_name)
        
        # Remove those sequences that are in the selected 46
        mAb_df = anti_merge_seqs(mAb_df, syn_seq_list)
        
        print("Starting to train models for "+mAb_name)

        # Train RFs
        RF_mAb_model, RF_mAb_test_data, RF_mAb_cv_score, RF_mAb_acc, RF_mAb_f1 = train_RF(mAb_df, parameters, 50, mAb_name, random_state, same_state = True)
        RF_cv_scores[i] = (RF_mAb_cv_score)
        RF_accuracy_scores[i] = (RF_mAb_acc)
        RF_f1_scores[i] = RF_mAb_f1

        # Train RNNs
        RNN_mAb_model, RNN_mAb_test_data, RNN_mAb_acc, RNN_mAb_f1 = train_RNN(mAb_df, mAb_name, random_state, same_state = True)
        RNN_accuracy_scores[i] = (RNN_mAb_acc)
        RNN_f1_scores[i] = (RNN_mAb_f1)
    
    RF_cv_scores.columns = mAb_names
    RF_accuracy_scores.columns = mAb_names
    RF_f1_scores.columns = mAb_names

    RNN_accuracy_scores.columns = mAb_names
    RNN_f1_scores.columns = mAb_names
    
    # save scores
    RF_cv_scores.to_csv('All_RF_Model_CVScores.tsv', sep = '\t', index = True)
    RF_accuracy_scores.to_csv('All_RF_Model_AccScores.tsv', sep = '\t', index = True)
    RF_f1_scores.to_csv('All_RF_Model_F1Scores.tsv', sep = '\t', index = True)
    
    RNN_accuracy_scores.to_csv('All_RNN_Model_AccScores.tsv', sep = '\t', index = True)
    RNN_f1_scores.to_csv('All_RNN_Model_F1Scores.tsv', sep = '\t', index = True)

"""
Run final function to train all models
"""
mAb_names = ['ACE2', 'LY16', 'LY555', 'REGN33', 'REGN87']
train_both_models(mAb_names, random_state)



























