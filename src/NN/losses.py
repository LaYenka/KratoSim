import numpy as np

# loss function and its derivative - express E = 1/N Sum(Y-Y*)^2, dE/dY = 2/N (Y-Y*)
def mse(y_true, y_pred):
    return np.mean(np.power(y_true-y_pred, 2));

def mse_prime(y_true, y_pred):
    return 2*(y_pred-y_true)/y_true.size;   
