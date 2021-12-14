# Import packages
import time
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import accuracy_score, matthews_corrcoef, f1_score, plot_precision_recall_curve, plot_confusion_matrix, plot_roc_curve
from sklearn.utils import shuffle
import joblib


"""
Additional code added by Beichen Gao for one-hot encoding, and matrix flattening
"""

def encode_onehot_padded(aa_seqs):
    '''
    one-hot encoding of a list of amino acid sequences with padding
    parameters:
        - aa_seqs : list with CDR3 sequences
    returns:
        - enc_aa_seq : list of np.ndarrays containing padded, encoded amino acid sequences
    '''
    ### Create an Amino Acid Dictionary
    aa_list = sorted(['A', 'C', 'D', 'E','F','G', 'H','I', 'K', 'L','M','N', \
              'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'])
    
    aa_dict = {char : l for char, l in zip(aa_list, np.eye(len(aa_list), k=0))}
    
    #####pad the longer sequences with '-' sign
    #1) identify the max length
    max_seq_len = max([len(x) for x in aa_seqs])
    #2) pad the shorter sequences with '-'
    aa_seqs = [seq + (max_seq_len - len(seq))*'-'
                      for i, seq in enumerate(aa_seqs)]
    
    # encode sequences:
    sequences=[]
    for seq in aa_seqs:
        e_seq=np.zeros((len(seq),len(aa_list)))
        count=0
        for aa in seq:
            if aa in aa_list:
                e_seq[count]=aa_dict[aa]
                count+=1
            else:
                print ("Unknown amino acid in peptides: "+ aa +", encoding aborted!\n")
        sequences.append(e_seq)
    enc_aa_seq = np.asarray(sequences)
    return enc_aa_seq


def flatten_matrix(encoded_aa):
    '''
    simple function to flatten 3D matrix of input AA list of dimensions 
    (data_size, positions, one-hot embedding) into 
    (data_size, flattened_embedding)
    '''
    mat = encoded_aa
    flat = np.reshape(mat, (mat.shape[0], -1)) #-1 infers last dimension
    return flat

def levenshtein_RBD(t, RBD):
    ''' 
    From Wikipedia article; Iterative with two matrix rows. 
    '''
    if RBD == t: 
        return 0
    elif len(RBD) == 0: 
        return len(t)
    elif len(t) == 0: 
        return len(RBD)
    v0 = [None] * (len(t) + 1)
    v1 = [None] * (len(t) + 1)
    for i in range(len(v0)):
        v0[i] = i
    for i in range(len(RBD)):
        v1[0] = i + 1
        for j in range(len(t)):
            cost = 0 if RBD[i] == t[j] else 1
            v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
        for j in range(len(v0)):
            v0[j] = v1[j]
            
    return v1[len(t)]


def make_dist_df(df):
    '''
    Takes input df of annotated sequences, and creates a dataframe of sequence counts by distance
    '''
    c = df.Distance.value_counts()
    c = pd.DataFrame(c).reset_index()
    c = c.rename(columns = {'index':'Distance', 'Distance':'Counts'})
    c.sort_values('Distance', inplace = True)
    return (c)

def clean_df(df, label, name, RBD_seq, threshold = 0, max_dist = 16):
    """
    Parameters
    ----------
    df : Pandas Dataframe
        Imported Pandas Dataframe read in from .tsv file with AA sequence in "junction_aa" column
        and consensus_count as the count column.
    threshold : Int
        Interger count threshold for removing noise (default is 5)
    max_dist : Int
        Maximum LD distance expected (default is 16)
    label : Int
        Label for sequence: 0 is negative, 1 is positive class
    name : String
        Name of the dataframe for saving processed data into a separate file and for the graph.

    Returns
    -------
    Cleaned dataframe, barplot of distances, and saves a copy of the cleaned dataframe as .tsv.
    
    Requires
    -------
    levenshtein_RBD function, make_dist_df function, matplotlib
    """
    # Filter dataframe by threshold
    a = df[df.consensus_count >= threshold]
    # Add Label to cleaned df
    a["Label"] = label
    
    # Add Distance to cleaned and original df for comparison
    a['Distance'] = a.apply(lambda row: levenshtein_RBD(row.junction_aa, RBD_seq), axis = 1)
    
    # Add 'Antibody' column
    a['Antibody'] = name
    
    # Keep only those with distance below max_dist threshold
    a = a[a.Distance <= max_dist]
    return (a)

def concat_dfs(df_list, name):
    """
    
    Parameters
    ----------
    df_list : List
        A list of dfs preprocessed using clean_df() to be combined for final dataset.
    name : String
        A string for the name of the final output tsv.

    Returns
    -------
    full concatenated dataframe, and saves it as a .tsv file.

    """
    
    final_df = pd.DataFrame()
    # Concatenate all dataframes in the list
    for df in df_list:
        final_df = pd.concat([final_df, df], axis = 0)
    
    # Remove duplicates
    final_df_cleaned = final_df.drop_duplicates(subset = 'junction_aa', keep = False)
    dropped_seqs = len(final_df) - len(final_df_cleaned)
    print('After cleaning',dropped_seqs,'Duplicated sequences were removed!')
    # Now separate by Label, and plot out distances
    
    minimum = int(min(final_df_cleaned.Distance))
    maximum = int(max(final_df_cleaned.Distance))
    graph_labels = list(map(str, range(minimum, (maximum+1))))
    graph_labels_df = pd.DataFrame(data = graph_labels, columns = {'Distance'})
    graph_labels_df['Distance'] = graph_labels_df['Distance'].astype(int) # Need to reset type to merge

    # Make positive and negative sets by label
    pos_set = make_dist_df(final_df_cleaned[final_df_cleaned.Label == 1])
    pos_set['Distance'] = pos_set['Distance'].astype(int)
    
    neg_set = make_dist_df(final_df_cleaned[final_df_cleaned.Label == 0])
    neg_set['Distance'] = neg_set['Distance'].astype(int)
    
    # Merge with graph_labels_df for final full dfs for graphing
    graph_pos = pd.merge(graph_labels_df, pos_set, how = 'left', on = 'Distance')
    graph_pos.fillna(0, inplace = True)
    
    graph_neg = pd.merge(graph_labels_df, neg_set, how = 'left', on = 'Distance')
    graph_neg.fillna(0, inplace = True)
    
    # Now prepare for graph
    pos_counts = graph_pos['Counts'].squeeze()
    neg_counts = graph_neg['Counts'].squeeze()
    
    a = np.arange(len(graph_labels))
    width = 0.35

    fig,ax = plt.subplots()
    
    rects1 = ax.bar(a-width/2, pos_counts, width, label = 'Positive', color = 'orange')
    rects2 = ax.bar(a+width/2, neg_counts, width, label = 'Negative', color = 'blue')

    ax.set_ylabel('Counts')
    ax.set_title(name+'Positive and Negative Sequence Counts by Distance')
    ax.set_xticks(a)
    ax.set_xticklabels(graph_labels)
    ax.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(name+'_combined_sequence_distribution')
    plt.close()

    # Finally, export cleaned dataframe and distance counts as a separate tsv
    final_df_cleaned.to_csv(name + '_combined.tsv', sep = '\t', index = False)
    graph_pos.to_csv(name + '_combined_poscounts.tsv', sep = '\t', index = False)
    graph_neg.to_csv(name + '_combined_negcounts.tsv', sep = '\t', index = False)
    
    return (final_df_cleaned)

def balance_edit_distance(df, distances, name, rand_state):
    dist = distances
    balanced_list = pd.DataFrame()
    for i in dist:
        pos = df[(df["Label"] == 1)]
        neg = df[(df["Label"] == 0)]
        a = pos.loc[pos['Distance'] == i]
        b = neg.loc[neg['Distance'] == i]
        if len(a) <= len(b):
            balanced_list = balanced_list.append(a)
            balanced_list = balanced_list.append(b.sample(len(a), random_state = rand_state))
        if len(b) < len(a):
            balanced_list = balanced_list.append(a.sample(len(b), random_state = rand_state))
            balanced_list = balanced_list.append(b)
            
    minimum = int(min(balanced_list.Distance))
    maximum = int(max(balanced_list.Distance))
    graph_labels = list(map(str, range(minimum, (maximum+1))))
    graph_labels_df = pd.DataFrame(data = graph_labels, columns = {'Distance'})
    graph_labels_df['Distance'] = graph_labels_df['Distance'].astype(int) # Need to reset type to merge

    # Make positive and negative sets by label
    pos_set = make_dist_df(balanced_list[balanced_list.Label == 1])
    pos_set['Distance'] = pos_set['Distance'].astype(int)
    
    neg_set = make_dist_df(balanced_list[balanced_list.Label == 0])
    neg_set['Distance'] = neg_set['Distance'].astype(int)
    
    # Merge with graph_labels_df for final full dfs for graphing
    graph_pos = pd.merge(graph_labels_df, pos_set, how = 'left', on = 'Distance')
    graph_pos.fillna(0, inplace = True)
    
    graph_neg = pd.merge(graph_labels_df, neg_set, how = 'left', on = 'Distance')
    graph_neg.fillna(0, inplace = True)
    
    # Now prepare for graph
    pos_counts = graph_pos['Counts'].squeeze()
    neg_counts = graph_neg['Counts'].squeeze()
    
    a = np.arange(len(graph_labels))
    width = 0.35

    fig,ax = plt.subplots()
    
    rects1 = ax.bar(a-width/2, pos_counts, width, label = 'Positive', color = 'orange')
    rects2 = ax.bar(a+width/2, neg_counts, width, label = 'Negative', color = 'blue')

    ax.set_ylabel('Counts')
    ax.set_title(name+' Balanced Positive and Negative Sequence Counts by Distance')
    ax.set_xticks(a)
    ax.set_xticklabels(graph_labels)
    ax.legend()

    fig.tight_layout()
    plt.show()
    plt.savefig(name+'_balanced_sequence_distribution')
    plt.close()
    
    # Add in a section where we keep the other sequences left out
    unseen_df = anti_merge(df, balanced_list)

    balanced_list.to_csv(name + '_balanced.tsv', sep = '\t', index = False)
    unseen_df.to_csv(name+'_unseen.tsv', sep = '\t', index = False)
    return(balanced_list)

def anti_merge(df_a, df_b):
    merge_df = pd.merge(df_a, df_b, how = 'outer', on = ['junction_aa', 'Label'], indicator = True)
    merge_df_left = merge_df[merge_df['_merge'] == 'left_only']
    merge_df_left.drop(['Distance_y', '_merge'], axis = 1, inplace = True)
    merge_df_left.rename(columns = {'Distance_x':'Distance'}, inplace = True)
    return (merge_df_left)

def one_hot(df):

    #extract length of longest sequence
    maxlen = max(df['junction_aa'].apply(len))

	#define alphabet (letters/nucleotides/AAs that occur)
    alphabet = 'ACDEFGHIKLMNPQRSTVWY*'
	#generate dictionary for label encoding
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))

	#prepare scaffold containing empty onehot matrices for each sequence
    input_data = np.zeros((len(df),maxlen,len(alphabet)))

	#for each sequence: encode each character as an integer
	# and based on that integer give the position a one hot vector of length #features (a,c,t,g)
    counter=0
    for sequence in df['junction_aa']:
		#get integer encoding (sequence --> integers)
        integer_encoded = [char_to_int[char] for char in sequence]
		#one hot encode integer encoded sequence
        onehot_encoded = list()
        for value in integer_encoded:
            letter = [0 for _ in range(len(alphabet))]
            letter[value] = 1
            onehot_encoded.append(letter)
		#padding (post)
		#pad each sequence that is shorter than maxlen
        while len(onehot_encoded) < maxlen:
            pad_array = [0 for _ in range(len(alphabet))]
            onehot_encoded.append(pad_array)
		#update scaffold entry
        input_data[counter]= onehot_encoded
        counter+=1		
		#print(sequence,onehot_encoded)
    return input_data

def RNN_encode_data(df):
    
    input_data = encode_onehot_padded(df.junction_aa)
    
    y = df.Label
    labels = y.to_numpy()
    
    return input_data, labels

def RNN_model_test_on_distances(list_of_distances, input_dataframe, model):
    """
    Function for testing model accuracy on one-hot encoded sequences at each distance

    Parameters
    ----------
    list_of_distances : List
        List of distances, for example [1, 2, 3, 4, 5].
    input_dataframe : DataFrame
        Dataframe of sequences, with aa sequences as a column titled 
        "junction_aa", class in "Labels" column, and distance in "Distance".
    model : Model
        Trained model, either RF or RNN should work.

    Returns
    -------
    output : DataFrame
        Outputs dataframe of scores for accuracy, and f1 by class and distance.

    """
    dist = list_of_distances
    data = input_dataframe
    # Create dataframe to append final scores into
    output = pd.DataFrame()
    pos_accuracies = []
    neg_accuracies = []
    total_accuracy = []
    shuff_accuracy = []
    total_f1 = []
    shuff_f1 = []
    distances = list(map(str, dist))
    for i in dist:
        # Extract relevant data from full dataframe
        data_dist = data[data.Distance == i]
        if len(data_dist) > 0:
            aa_seqs, y_dist = RNN_encode_data(data_dist)
            y_shuff = shuffle(y_dist)
            # Use model and predict
            pred_dist = model.predict(aa_seqs)
            pred_dist[pred_dist <= 0.5] = 0
            pred_dist[pred_dist > 0.5] = 1
            # Calculate scores
            score_dist = accuracy_score(pred_dist, y_dist)
            score_dist_shuff = accuracy_score(pred_dist, y_shuff)
            f1_dist = f1_score(pred_dist, y_dist)
            f1_dist_shuff = f1_score(pred_dist, y_shuff)
            # Append to dfs
            total_accuracy.append(score_dist)
            shuff_accuracy.append(score_dist_shuff)
            total_f1.append(f1_dist)
            shuff_f1.append(f1_dist_shuff)
        elif len(data_dist) == 0:
            score_dist = 0
            score_dist_shuff = 0
            f1_dist = 0
            f1_dist_shuff = 0
            total_accuracy.append(score_dist)
            shuff_accuracy.append(score_dist_shuff)
            total_f1.append(f1_dist)
            shuff_f1.append(f1_dist_shuff)
        # Now for POSITIVES
        data_dist_pos = data_dist[(data_dist['Label'] == 1)]
        if len(data_dist_pos) > 0:
            aa_seqs_pos, y_dist_pos = RNN_encode_data(data_dist_pos)
            # Use model and predict 
            pred_dist_pos = model.predict(aa_seqs_pos)
            pred_dist_pos[pred_dist_pos <= 0.5] = 0
            pred_dist_pos[pred_dist_pos > 0.5] = 1
            # Calculate scores
            score_pos = accuracy_score(pred_dist_pos,y_dist_pos)
            # Append to dataframe
            pos_accuracies.append(score_pos)
        elif len(data_dist_pos) == 0:
            score_pos = 0
            pos_accuracies.append(score_pos)
        # Now for NEGATIVES
        data_dist_neg = data_dist[(data_dist['Label'] == 0)]
        if len(data_dist_neg) > 0:
            aa_seqs_neg, y_dist_neg = RNN_encode_data(data_dist_neg)
            # Use model to predict
            pred_dist_neg = model.predict(aa_seqs_neg)
            pred_dist_neg[pred_dist_neg <= 0.5] = 0
            pred_dist_neg[pred_dist_neg > 0.5] = 1
            # Calculate Scores
            score_neg = accuracy_score(pred_dist_neg, y_dist_neg)
            # Append to dataframe
            neg_accuracies.append(score_neg)
        elif len(data_dist_neg) == 0:
            score_neg = 0
            neg_accuracies.append(score_neg)
    output['Pos_Accuracies'] = pos_accuracies
    output['Neg_Accuracies'] = neg_accuracies
    output['Total_Acc'] = total_accuracy
    output['Shuff_Acc'] = shuff_accuracy
    output['Total_F1'] = total_f1
    output['Shuff_F1'] = shuff_f1
    output.index = distances
    return (output)

def RF_model_test_on_distances(list_of_distances, input_dataframe, model):
    dist = list_of_distances
    data = input_dataframe
    # Create dataframe to append final scores into
    output = pd.DataFrame()
    pos_accuracies = []
    neg_accuracies = []
    total_accuracy = []
    shuff_accuracy = []
    total_f1 = []
    shuff_f1 = []
    distances = list(map(str, dist))
    for i in dist:
        # Extract relevant data from full dataframe
        data_dist = data[data.Distance == i]
        if len(data_dist) > 0:
            aa_seqs = flatten_matrix(encode_onehot_padded(data_dist.junction_aa))
            y_dist = data_dist.Label
            y_shuff = shuffle(y_dist)
            # Use model and predict
            pred_dist = model.predict(aa_seqs)
            # Calculate scores
            score_dist = accuracy_score(pred_dist, y_dist)
            score_dist_shuff = accuracy_score(pred_dist, y_shuff)
            f1_dist = f1_score(pred_dist, y_dist)
            f1_dist_shuff = f1_score(pred_dist, y_shuff)
            # Append to dfs
            total_accuracy.append(score_dist)
            shuff_accuracy.append(score_dist_shuff)
            total_f1.append(f1_dist)
            shuff_f1.append(f1_dist_shuff)
        elif len(data_dist) == 0:
            score_dist = 0
            score_dist_shuff = 0
            f1_dist = 0
            f1_dist_shuff = 0
            total_accuracy.append(score_dist)
            shuff_accuracy.append(score_dist_shuff)
            total_f1.append(f1_dist)
            shuff_f1.append(f1_dist_shuff)
        # Now for POSITIVES
        data_dist_pos = data_dist[(data_dist['Label'] == 1)]
        if len(data_dist_pos) > 0:
            aa_seqs_pos = flatten_matrix(encode_onehot_padded(data_dist_pos.junction_aa))
            y_dist_pos = data_dist_pos.Label
            # Use model and predict 
            pred_dist_pos = model.predict(aa_seqs_pos)
            # Calculate scores
            score_pos = accuracy_score(pred_dist_pos,y_dist_pos)
            # Append to dataframe
            pos_accuracies.append(score_pos)
        elif len(data_dist_pos) == 0:
            score_pos = 0
            pos_accuracies.append(score_pos)
        # Now for NEGATIVES
        data_dist_neg = data_dist[(data_dist['Label'] == 0)]
        if len(data_dist_neg) > 0:
            aa_seqs_neg = flatten_matrix(encode_onehot_padded(data_dist_neg.junction_aa))
            y_dist_neg = data_dist_neg.Label
            # Use model to predict
            pred_dist_neg = model.predict(aa_seqs_neg)
            # Calculate Scores
            score_neg = accuracy_score(pred_dist_neg, y_dist_neg)
            # Append to dataframe
            neg_accuracies.append(score_neg)
        elif len(data_dist_neg) == 0:
            score_neg = 0
            neg_accuracies.append(score_neg)
    output['Pos_Accuracies'] = pos_accuracies
    output['Neg_Accuracies'] = neg_accuracies
    output['Total_Acc'] = total_accuracy
    output['Shuff_Acc'] = shuff_accuracy
    output['Total_F1'] = total_f1
    output['Shuff_F1'] = shuff_f1
    output.index = distances
    return (output)





