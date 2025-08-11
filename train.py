# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 12:49:13 2025

@author: ctbac
"""

import numpy as np
from sklearn.model_selection import train_test_split
from tensorflow.keras.utils import to_categorical
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Dropout, Flatten, Dense, InputLayer
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ReduceLROnPlateau
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects import r
import numpy as np
from rpy2.robjects import numpy2ri
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from tensorflow.keras.utils import plot_model
from collections import Counter
from sklearn.metrics import precision_score, recall_score
import random



def expand_by_tens(x):
    result = []
    for xi in x:
        start = int(xi * 10)
        end = int((xi + 1) * 10)
        result.extend(range(start, end))
    return np.array(result)

def count_elements_in_bins(arr, bin_size=5):
    arr = np.asarray(arr)
    bins = (arr // bin_size) * bin_size
    count = Counter(bins)
    result = {f"[{k}-{k + bin_size - 1}]": v for k, v in sorted(count.items())}
    return result

def count_bins_exceeding_threshold(binned_counts: dict, c: int) -> int:
    return sum(1 for count in binned_counts.values() if count > c)


def most_common_in_chunks(x, c):
    x = np.asarray(x)
    n = len(x)
    y = []

    for i in range(0, n, c):
        chunk = x[i:i + c]
        if len(chunk) == 0:
            continue
        counter = Counter(chunk)
        most_common_element = counter.most_common(1)[0][0]
        y.append(most_common_element)
    
    return np.array(y)

r['load']('res.RData')
res = robjects.globalenv['res']
kep11 = np.array(res[0])
kep15 = np.array(res[1])
kep19 = np.array(res[2])
kep11smo = np.array(res[3])
kep15smo = np.array(res[4])
kep19smo = np.array(res[5])
pg = np.array(res[6])
pgsmo = np.array(res[7])
#ABCDE

acc = []
pre = []
rcl = []
random.seed(1)
for ep in range(10):
    # ind = np.concatenate([np.arange(0, 1000),np.arange(3000, 4000), np.arange(4000, 5000)]) # A-D-E
    # y = np.array([0]*1000 + [1]*1000 + [2]*1000) 
    
    ind = np.concatenate([np.arange(0, 1000), np.arange(4000, 5000)]) # A-E
    y = np.array([0]*1000 + [1]*1000 )  

    # ind = np.concatenate([np.arange(0, 4000), np.arange(4000, 5000)]) # ABCD-E
    # y = np.array([0]*4000 + [1]*1000 )  
    
    # ind = np.concatenate([np.arange(0, 4000), np.arange(4000, 5000)]) #  AB-CD-E
    # y = np.array([0]*2000 + [1]*2000 + [2]*1000) 

    X = kep11[ind,]
    
    # Step 1: data
    X = np.expand_dims(X, axis=-1)               # (5000, 99, 46, 1)
    y_cat = to_categorical(y, num_classes=2)     # (5000, 4)
    
    indices = np.arange(X.shape[0]/10)
    train_idx, test_idx = train_test_split(indices, test_size=0.25, random_state=1, shuffle=True)#stratify=y_cat[::10])
    train_idx = expand_by_tens(train_idx)
    test_idx = expand_by_tens(test_idx)
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y_cat[train_idx], y_cat[test_idx]
    
    #X_train, X_test, y_train, y_test = train_test_split(X, y_cat, test_size=0.2, stratify=y)
    
    # Step 2: model
    model = Sequential()
    model.add(InputLayer(input_shape=(99, 46, 1)))
    
    model.add(Conv2D(filters=32, kernel_size=(3, 3), activation='relu', padding='same'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    
    model.add(Conv2D(filters=64, kernel_size=(3, 3), activation='relu', padding='same'))
    model.add(MaxPooling2D(pool_size=(2, 2)))
    model.add(Dropout(0.25))
    
    
    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.25))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(2, activation='softmax'))
    
    # Step 3
    optimizer = Adam(learning_rate=1e-3)
    model.compile(optimizer=optimizer,
                  loss='categorical_crossentropy',
                  metrics=['accuracy'])
    
    # Step 4: (Reduce on Plateau)
    lr_scheduler = ReduceLROnPlateau(
        monitor='val_loss',
        factor=0.5,           # 学习率乘以 0.5
        patience=3,           # 如果验证集3轮没提升，则衰减
        min_lr=5e-6,          # 最小学习率
        verbose=1
    )
    
    # Step 5: fit
    history = model.fit(
        X_train, y_train,
        batch_size=32,
        epochs=40,
        validation_data=(X_test, y_test),
        callbacks=[lr_scheduler]
    )
    
    # Step 6: plot curve
    # plt.figure(figsize=(8, 5))
    # plt.plot(history.history['accuracy'], label='Training Accuracy', marker='o')
    # plt.plot(history.history['val_accuracy'], label='Validation Accuracy', marker='s')
    # plt.title('Accuracy Curve')
    # plt.xlabel('Epoch')
    # plt.ylabel('Accuracy')
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    
    
    y_true = np.argmax(y_test, axis=1)
    y_pred = model.predict([X_test])
    y_pred_classes = np.argmax(y_pred, axis=1)
    
    # confusion matrix
    y_true = most_common_in_chunks(y_true, 10)
    y_pred_classes = most_common_in_chunks(y_pred_classes, 10)
    cm = confusion_matrix(y_true, y_pred_classes)
    # disp = ConfusionMatrixDisplay(confusion_matrix=cm)
    # disp.plot(cmap='Blues', values_format='d')  # 格式化为整数
    # plt.title("Confusion Matrix")
    # plt.show()
    wrong_ind = np.where(y_true!=y_pred_classes)[0]
    # print('Number of misclassification: ', len(wrong_ind))
    precision = precision_score(y_true, y_pred_classes, average='macro')
    recall = recall_score(y_true, y_pred_classes, average='macro')
    
    acc.append(1-len(wrong_ind)/len(y_true))
    pre.append(precision)
    rcl.append(recall)

print(np.mean(acc))
print(np.mean(pre))
print(np.mean(rcl))
print(np.std(acc))
print(np.std(pre))
print(np.std(rcl))

# import shutil
# import os
# import sys
# graphviz_path = r"C:\Program Files\Graphviz\bin"
# os.environ["PATH"] += os.pathsep + graphviz_path
# print("dot in PATH (after manual add):", shutil.which("dot"))

# plot_model(model, to_file='model_structure.png', show_shapes=True, show_layer_names=True)
# model.save("cnn_model.h5")
