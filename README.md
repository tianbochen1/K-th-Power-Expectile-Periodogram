# K-th-Power-Expectile-Periodogram

This paper introduces the $k$-th power expectile periodogram (KEP) as a generalization of the quantile periodogram (QP) and expectile periodogram (EP) for time series analysis. Constructed from trigonometric $k$-th power expectile regression (KER), the KEP provides a flexible tool to balance robustness and efficiency. The KEP retains the key properties of the ordinary periodogram as a frequency-domain representation of serial dependence in time series, while offering a more comprehensive understanding by examining the data across the entire range of expectile levels with power $1<k\leq2$. Simulations demonstrate the capability of the KEP in detecting hidden periodicities, while maintaining robustness against outliers. We establish the asymptotic theory and investigate the relationship between the KEP and the so-called $k$-th power expectile spectrum (KES). Finally, we leverage the inherent two-dimensional property of the KEP to train a convolutional neural network (CNN) to classify the epilepsy electroencephalogram (EEG) data, where the proposed estimator outperforms competitive estimators.

One can download the code and data to reproduce the results in the paper uploaded.

Code:
"fn.r" contains functions to compute the KEP.
"eeg example.r" and "sp example.r" contains the code to reprodeuce the figures in Section 2.3.
"kepsim.r" reproduces the results in the simulations.


Data:

"Patient-3_seizure-2_Q5-Q6_17-21.mat", EEG data.
"sp.csv", S&P500 data.

For further inquiries, please contact Tianbo Chen (chentianbo@ahu.edu.cn).
