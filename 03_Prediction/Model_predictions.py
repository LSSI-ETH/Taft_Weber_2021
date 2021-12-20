#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:54:27 2021

@author: bgao

Code for "Predictive profiling of SARS-CoV-2 variants by deep mutational learning"
Taft J and Weber C.

Using trained RF (after calibration - update names to select best ones) and RNN models
to predict on new prospective sequences generated from current VoC
"""

"""
Load all functions needed
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
Importing all models for future predictions
- Note, in our case, the isotonic models were the best performing after calibration

"""
mAb_names_list = ['ACE2', 'LY16', 'LY555', 'REGN33', 'REGN87']
RF_models_list = []
RNN_models_list = []

def mAb_RF_models_import(mAb_name):
    # Import models
    RNN = keras.models.load_model(mAb_name+'_RNN_All_3LSTM.h5')
    RF = joblib.load(mAb_name+'_isotonic_RFmodel.gz')
    RNN_models_list.append(RNN)
    RF_models_list.append(RF)
    
for mAb_name in mAb_names_list:
    mAb_RF_models_import(mAb_name)
    
RNN_labels = ['RNN_ACE2', 'RNN_LY16', 'RNN_LY555', 'RNN_REGN33', 'RNN_REGN87']

RF_labels = ['RF_ACE2', 'RF_LY16', 'RF_LY555', 'RF_REGN33', 'RF_REGN87']

"""
Now predict on list of variant sequences from the same directory
with naming scheme: all_seqs_dist123_valid_(variant_name).tsv

TSV should contain column of variant aa sequences, label for distance from VoC, and
positions of mutations
"""
variants = ['alpha', 'beta', 'gamma', 'kappa', 'wuhan', 'new_variant']

for variant_name in variants:
    variant_df = pd.read_table('all_seqs_dist123_valid_'+variant_name+".tsv")
    X_RF = flatten_matrix(encode_onehot_padded(variant_df.sequence_aa))
    X_RNN = encode_onehot_padded(variant_df.sequence_aa)
    # Predict using RNN
    for i in range(len(RNN_labels)):
        col_label = RNN_labels[i]
        model = RNN_models_list[i]
        variant_df[col_label] = model.predict(X_RNN)
    # Predict using RF
    for i in range(len(RF_labels)):
        col_label = RF_labels[i]
        model = RF_models_list[i]
        variant_df[col_label] = model.predict_proba(X_RF)[:,1]
    # Save final df with predictions
    variant_df.to_csv(variant_name+"model_predictions.csv")

