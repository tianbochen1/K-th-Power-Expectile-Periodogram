# K-th-Power-Expectile-Periodogram

This work proposes the $k$-th power expectile periodogram (KEP) as a generalization of the quantile periodogram (QP) and expectile periodogram (EP) for time series analysis. By raising the check function to the $k$-th power, the KEP provides a flexible tool to balance robustness and efficiency. The KEP offers a more comprehensive understanding by examining the time series across the entire range of expectile levels at power $1<k\leq2$. Simulations demonstrate the capability of the proposed periodogram in detecting hidden periodicities, while maintaining robustness against outliers. We establish the asymptotic theory and investigate the relationship between the estimator and the so-called $k$-th power expectile spectrum (KES). Finally, we leverage the inherent two-dimensional property of the KEP, and train a convolution neural network (CNN) to classify the epileptic electroencephalogram (EEG) data, where the proposed estimator outperforms the ordinary periodogram.

One can download the code and data to reproduce the results in the paper uploaded.

Code:
"fn.r" contains functions to compute the KEP.
"eeg example.r" and "sp example.r" contains the code to reprodeuce the figures in Section 2.3.
"kepsim.r" reproduces the results in the simulations.
"train.py" train the CNN model in Section 5.

Data:

"Patient-3_seizure-2_Q5-Q6_17-21.mat", EEG data.
"sp.csv", S&P500 data.

For further inquiries, please contact Tianbo Chen (chentianbo@ahu.edu.cn).
