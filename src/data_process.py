
import numpy as np
import pickle
from sklearn.model_selection import KFold
import torch
import pandas as pd
import numpy as np
from sklearn import preprocessing


def csv_load(filename, cost_from_file):
    # csv format
    # 1st row : cost from 2nd col (align with column name in 2nd row) if
    # cost_from_file if True. Else below description of 2nd row is for 1st row
    # 2nd row : columns name starting with 'label' followed by features' name
    # 3rd ~ rows : label and feature values
    header = 1 if cost_from_file else 0
    df = pd.read_csv(filename, header=header)
    labels = df['label'].values.astype(np.int)
    df = df.iloc[:, 2:] # assume 1st and 2nd col are id,label
    #exist = np.where(pd.isnull(df), 1, 0).astype(np.uint8)
    exist = np.where(pd.isnull(df), 0, 1).astype(np.uint8)

    def norm_to_zero_one(df):
        return (df - df.min()) * 1.0 / (df.max() - df.min())
    def std_norm(df):
        return (df - df.mean()) / df.std()

    df = norm_to_zero_one(df)
    df = df.fillna(0)#method='backfill')

    if cost_from_file:
        cost = pd.read_csv(filename, nrows=1,header=None)
        cost = cost.values.reshape(-1)[2:]
        assert len(cost) == df.shape[1]
    else:
        cost = None

    return df.values.astype(np.float32), exist, labels, cost

def data_split_n_ready(features, exist, labels,
        mode='cv', random_seed=123,
        val_test_split=np.array([0.25, 0.25]), action2features=None,
        shuffle=True):
    dataset_size = len(features)
    indices = list(range(dataset_size))
    assert np.sum(val_test_split) < 1
    split = np.floor(val_test_split * dataset_size).astype(np.int)
    split = np.cumsum(split)
    if shuffle:
        np.random.seed(random_seed)
        np.random.shuffle(indices)

    train_indices = indices[split[1]:]
    val_indices = indices[:split[0]]
    test_indices = indices[split[0]:split[1]]

    def pick(indices, iter):
        exist_ = exist[indices] if exist is not None else None
        return DataTemp(features[indices], labels[indices], exist_,
                n_classes=(np.amax(labels) + 1), # label from 0 to n_classes-1
                action2features=action2features, iter=iter)

    return pick(train_indices, True), \
            pick(val_indices, False), pick(test_indices, False)

def data_load(data_type='csv', # csv
              random_seed=123,
              n_datapoints=20000, # ignored when data_type is csv
              csv_filename=None,
              action2features=None,
              val_test_split=np.array([0.25, 0.25]),
              cost_from_file=False):
    if data_type == 'csv':
        features, exist, labels, cost = csv_load(csv_filename, cost_from_file)

    if cost_from_file == False:
        n_features = 130
        cost = None
        cost_uniform = -abs(np.random.uniform(low=-50,high=-1,size=n_features))
        #cost_normal = -abs(np.random.normal(loc=-30,scale=10,size=n_features))
        #cost_exp = -abs(np.random.exponential(scale=10,size=n_features))
        cost = cost_uniform #cost_normal  cost_exp 

    a,b,c = data_split_n_ready(features, exist, labels,
                            val_test_split=val_test_split,
                            action2features=action2features)
    return a,b,c,cost



#### cost simulation
#    cost = None
#    cost_uniform = -abs(np.random.uniform(low=-50,high=-1,size=n_features))
#    cost_normal = -abs(np.random.normal(loc=-30,scale=10,size=n_features))
#    cost_exp = -abs(np.random.exponential(scale=10,size=n_features))
#    cost = cost_uniform

class DataTemp:
    def __init__(self, features, labels, exist, n_classes, shuffle=True, iter=True,
            action2features=None):
        self.features = features
        self.labels = labels
        self.exist = exist
        self.shuffle = shuffle
        self.n_data, self.n_features = features.shape
        self.n_classes = n_classes
        self.index = 0
        self.iter = iter
        self.action2features = action2features
        self.n_actions = len(action2features) + 1 if action2features is not None \
                else features.shape[1] + 1

    def next_batch(self, batch_size):
        if iter:
            new_index = (self.index + batch_size) % self.n_data
        else:
            if self.index == self.n_data:
                return None # Done
            new_index = min(self.index + batch_size, self.n_data)

        if self.index + batch_size <= self.n_data:
            features = self.features[self.index: self.index + batch_size]
            labels = self.labels[self.index: self.index + batch_size]
            exist = self.exist[self.index: self.index + batch_size] \
                if self.exist is not None else None
        else:
            features = self.features[self.index:]
            labels = self.labels[self.index:]
            exist = self.exist[self.index:] if self.exist is not None else None
            if self.iter:
                if self.shuffle:
                    p = np.random.permutation(self.n_data)
                    self.features = self.features[p]
                    self.labels = self.labels[p]
                    self.exist = self.exist[p] if self.exist is not None else None
                features = np.concatenate((features,
                    self.features[:new_index]), axis=0)
                labels = np.concatenate((labels,
                    self.labels[:new_index]), axis=0)
                exist = np.concatenate((exist, self.exist[:new_index]), axis=0) \
                        if self.exist is not None else None
        self.index = new_index
        return features, labels, exist
