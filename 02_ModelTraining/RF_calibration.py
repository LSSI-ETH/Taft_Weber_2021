#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:26:21 2021

@author: bgao

Code for "Predictive profiling of SARS-CoV-2 variants by deep mutational learning"
Taft J and Weber C.

Calibrate RF models trained on ACE2/mAb binding and escape sequences
"""

import os
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, matthews_corrcoef, f1_score
from sklearn.utils import shuffle
import h5py

import seaborn as sns
import matplotlib.pyplot as plt

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

import joblib

from helpers import read_blosum, encode_blosum, encode_onehot_padded, flatten_matrix, levenshtein_RBD, make_dist_df, clean_df, concat_dfs, balance_edit_distance, anti_merge, make_RF_model, model_test_on_distances

"""
Create function to evaluate calibrated models by distance
"""

from sklearn.calibration import CalibratedClassifierCV

def eval_by_dist(distances, model, test_df, mAb_name, model_name):
    model_performance = model_test_on_distances(distances, test_df, model)

    graph_labels = list(map(str, distances))
    graph_positive = model_performance['Pos_Accuracies'].squeeze()
    graph_negative = model_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,ax = plt.subplots()
    rects1 = ax.bar(a-width/2, graph_positive, width, label = 'Positive', color = 'orange')
    rects2 = ax.bar(a+width/2, graph_negative, width, label = 'Negative', color = 'blue')

    ax.set_ylabel('Accuracies')
    ax.set_title(mAb_name + model_name + ' Model Accuracy on Test Sequences by Edit Distance')
    ax.set_xticks(a)
    ax.set_xticklabels(graph_labels)
    ax.legend()
    fig.tight_layout()
    plt.show()
    plt.savefig(mAb_name + model_name + '_Model_Accuracy_by_Distance')
    plt.close()
    model_performance.to_csv(mAb_name + model_name + 'RFModel_Accuracies_by_Distance.csv')
    
    # Then predict accuracies on unseen data
    unseen_data = pd.read_table(mAb_name+'_RF_unseen.tsv')
    unseen_performance = model_test_on_distances(distances, unseen_data, model)

    graph_labels = list(map(str, distances))
    graph_positive2 = unseen_performance['Pos_Accuracies'].squeeze()
    graph_negative2 = unseen_performance['Neg_Accuracies'].squeeze()

    a = np.arange(len(graph_labels))
    width = 0.35

    fig,dx = plt.subplots()
    rects1 = dx.bar(a-width/2, graph_positive2, width, label = 'Positive', color = 'orange')
    rects2 = dx.bar(a+width/2, graph_negative2, width, label = 'Negative', color = 'blue')

    dx.set_ylabel('Accuracies')
    dx.set_title(mAb_name + model_name + ' Model Unseen Accuracy on Test Sequences by Edit Distance')
    dx.set_xticks(a)
    dx.set_xticklabels(graph_labels)
    dx.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(mAb_name + model_name + '_Model_Unseen_Accuracy_by_Distance')
    plt.close()
    unseen_performance.to_csv(mAb_name+model_name+'RFModel_Unseen_Accuracies_by_Distance.csv')
    
    #Finally, make the PR-ROC curve
    from sklearn.metrics import roc_curve, auc
    
    X_test_oh = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    y_test_shuff = shuffle(y_test)
    
    y_pred_RF = model.predict_proba(X_test_oh)[:,1].ravel()
    fpr_RF, tpr_RF, thresholds_RF = roc_curve(y_test, y_pred_RF)
    auc_RF = auc(fpr_RF, tpr_RF)
    
    fpr_RF_shuff, tpr_RF_shuff, thresholds_RF_shuff = roc_curve(y_test_shuff, y_pred_RF)
    auc_RF_shuff = auc(fpr_RF_shuff, tpr_RF_shuff)
    
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(fpr_RF, tpr_RF, label=(mAb_name+model_name+'_RF (Full) (area = {:.3f})').format(auc_RF))
    plt.plot(fpr_RF_shuff, tpr_RF_shuff, label=(mAb_name+model_name+'_RF (Shuffled) (area = {:.3f})').format(auc_RF_shuff))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve for Full Models')
    plt.legend(loc='best')
    plt.show()
    plt.savefig(mAb_name+model_name+"_RF_ROCCurve")
    plt.close()

"""
Create functions to train calibrated models
Note: Currently this is done in a roundabout way, models are retrained using
best parameters found from RandomSearchCV during training, and then calibrated using
CalibratedClassifierCV

Best parameters are saved in the "mAb_RF_GridSearchCV.csv" for each mAb

Future versions will have this calibration wrapped in the model training step.
"""

def train_ACE2_calibrated(mAb_name):
    # Load in data
    mAb_model = joblib.load(mAb_name+"_RF_model.gz")
    train_df = pd.read_csv(mAb_name+"_RF_train_data.csv")
    test_df = pd.read_csv(mAb_name+"_RF_test_data.csv")
    
    # Make x_train and y_train
    X_train = flatten_matrix(encode_onehot_padded(train_df.junction_aa))
    y_train = train_df.Label
    X_test = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    
    # Make calibrated models
    RF_model = RandomForestClassifier(n_estimators = 400, max_features = 'sqrt', max_depth = 60)
    clf_isotonic = CalibratedClassifierCV(RF_model, cv = 5, method = 'isotonic')
    clf_sigmoid = CalibratedClassifierCV(RF_model, cv = 5, method = 'sigmoid')

    # Fit models
    clf_isotonic.fit(X_train, y_train)
    clf_sigmoid.fit(X_train, y_train)
    
    # Evaluate models
    distances = list(range(1,11))
    
    eval_by_dist(distances, clf_isotonic, test_df, mAb_name, "_isotonic")
    eval_by_dist(distances, clf_sigmoid, test_df, mAb_name, "_sigmoid")
    
    # Save models
    joblib.dump(clf_isotonic, mAb_name+"_isotonic_RFmodel.gz")
    joblib.dump(clf_sigmoid, mAb_name+"_sigmoid_RFmodel.gz")
    return mAb_model, clf_isotonic, clf_sigmoid

def train_LY16_calibrated(mAb_name):
    # Load in data
    mAb_model = joblib.load(mAb_name+"_RF_model.gz")
    train_df = pd.read_csv(mAb_name+"_RF_train_data.csv")
    test_df = pd.read_csv(mAb_name+"_RF_test_data.csv")
    
    # Make x_train and y_train
    X_train = flatten_matrix(encode_onehot_padded(train_df.junction_aa))
    y_train = train_df.Label
    X_test = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    
    # Make calibrated models
    RF_model = RandomForestClassifier(n_estimators = 450, max_features = 'sqrt', max_depth = 150)
    clf_isotonic = CalibratedClassifierCV(RF_model, cv = 5, method = 'isotonic')
    clf_sigmoid = CalibratedClassifierCV(RF_model, cv = 5, method = 'sigmoid')

    # Fit models
    clf_isotonic.fit(X_train, y_train)
    clf_sigmoid.fit(X_train, y_train)
    
    # Evaluate models
    distances = list(range(1,11))
    
    eval_by_dist(distances, clf_isotonic, test_df, mAb_name, "_isotonic")
    eval_by_dist(distances, clf_sigmoid, test_df, mAb_name, "_sigmoid")
    
    # Save models
    joblib.dump(clf_isotonic, mAb_name+"_isotonic_RFmodel.gz")
    joblib.dump(clf_sigmoid, mAb_name+"_sigmoid_RFmodel.gz")
    return mAb_model, clf_isotonic, clf_sigmoid

def train_LY555_calibrated(mAb_name):
    # Load in data
    mAb_model = joblib.load(mAb_name+"_RF_model.gz")
    train_df = pd.read_csv(mAb_name+"_RF_train_data.csv")
    test_df = pd.read_csv(mAb_name+"_RF_test_data.csv")
    
    # Make x_train and y_train
    X_train = flatten_matrix(encode_onehot_padded(train_df.junction_aa))
    y_train = train_df.Label
    X_test = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    
    # Make calibrated models
    RF_model = RandomForestClassifier(n_estimators = 300, max_features = 'auto', max_depth = 80)
    clf_isotonic = CalibratedClassifierCV(RF_model, cv = 5, method = 'isotonic')
    clf_sigmoid = CalibratedClassifierCV(RF_model, cv = 5, method = 'sigmoid')

    # Fit models
    clf_isotonic.fit(X_train, y_train)
    clf_sigmoid.fit(X_train, y_train)
    
    # Evaluate models
    distances = list(range(1,11))
    
    eval_by_dist(distances, clf_isotonic, test_df, mAb_name, "_isotonic")
    eval_by_dist(distances, clf_sigmoid, test_df, mAb_name, "_sigmoid")
    
    # Save models
    joblib.dump(clf_isotonic, mAb_name+"_isotonic_RFmodel.gz")
    joblib.dump(clf_sigmoid, mAb_name+"_sigmoid_RFmodel.gz")
    return mAb_model, clf_isotonic, clf_sigmoid

def train_REGN33_calibrated(mAb_name):
    # Load in data
    mAb_model = joblib.load(mAb_name+"_RF_model.gz")
    train_df = pd.read_csv(mAb_name+"_RF_train_data.csv")
    test_df = pd.read_csv(mAb_name+"_RF_test_data.csv")
    
    # Make x_train and y_train
    X_train = flatten_matrix(encode_onehot_padded(train_df.junction_aa))
    y_train = train_df.Label
    X_test = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    
    # Make calibrated models
    RF_model = RandomForestClassifier(n_estimators = 400, max_features = 'auto', max_depth = 80)
    clf_isotonic = CalibratedClassifierCV(RF_model, cv = 5, method = 'isotonic')
    clf_sigmoid = CalibratedClassifierCV(RF_model, cv = 5, method = 'sigmoid')

    # Fit models
    clf_isotonic.fit(X_train, y_train)
    clf_sigmoid.fit(X_train, y_train)
    
    # Evaluate models
    distances = list(range(1,11))
    
    eval_by_dist(distances, clf_isotonic, test_df, mAb_name, "_isotonic")
    eval_by_dist(distances, clf_sigmoid, test_df, mAb_name, "_sigmoid")
    
    # Save models
    joblib.dump(clf_isotonic, mAb_name+"_isotonic_RFmodel.gz")
    joblib.dump(clf_sigmoid, mAb_name+"_sigmoid_RFmodel.gz")
    return mAb_model, clf_isotonic, clf_sigmoid

def train_REGN87_calibrated(mAb_name):
    # Load in data
    mAb_model = joblib.load(mAb_name+"_RF_model.gz")
    train_df = pd.read_csv(mAb_name+"_RF_train_data.csv")
    test_df = pd.read_csv(mAb_name+"_RF_test_data.csv")
    
    # Make x_train and y_train
    X_train = flatten_matrix(encode_onehot_padded(train_df.junction_aa))
    y_train = train_df.Label
    X_test = flatten_matrix(encode_onehot_padded(test_df.junction_aa))
    y_test = test_df.Label
    
    # Make calibrated models
    RF_model = RandomForestClassifier(n_estimators = 150, max_features = 'sqrt', max_depth = None)
    clf_isotonic = CalibratedClassifierCV(RF_model, cv = 5, method = 'isotonic')
    clf_sigmoid = CalibratedClassifierCV(RF_model, cv = 5, method = 'sigmoid')

    # Fit models
    clf_isotonic.fit(X_train, y_train)
    clf_sigmoid.fit(X_train, y_train)
    
    # Evaluate models
    distances = list(range(1,11))
    
    eval_by_dist(distances, clf_isotonic, test_df, mAb_name, "_isotonic")
    eval_by_dist(distances, clf_sigmoid, test_df, mAb_name, "_sigmoid")
    
    # Save models
    joblib.dump(clf_isotonic, mAb_name+"_isotonic_RFmodel.gz")
    joblib.dump(clf_sigmoid, mAb_name+"_sigmoid_RFmodel.gz")
    return mAb_model, clf_isotonic, clf_sigmoid

def train_all_calibrated_RF():
    train_ACE2_calibrated('ACE2')
    train_LY16_calibrated('LY16')
    train_LY555_calibrated('LY555')
    train_REGN33_calibrated('REGN33')
    train_REGN87_calibrated('REGN87')

train_all_calibrated_RF()


