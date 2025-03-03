# -*- coding: utf-8 -*-
import numpy as np

def gradient(params, y, x1):
    n = len(y)
    y = np.array(y)
    x1 = np.array(x1)
    a0 = params[0]
    b1 = params[1]
    grad = np.array([
    -2*(y - a0 - b1*x1).sum()/n,
    -2*(np.dot(np.transpose(y - a0 - b1*x1), x1)).sum()/n
    ])
    return grad
