import numpy as np
import matplotlib.pyplot as plt

def con_grad(A,x,b):
    r = b - A @ x
    p = r
    j = 0

    conv = []
    err = []

    conv.append(np.linalg.norm(r))

    if b.all() == 0:
        err.append(np.linalg.norm(x))

    while np.linalg.norm(r)/np.linalg.norm(b) > 1e-14 and j < 500:
        alpha = (r.T @ r) / (p.T @ A @ p)
        x = x + alpha * p
        rnew = r - alpha * A @ p
        beta = (rnew.T @ rnew)/(r.T @ r)
        p = rnew + beta * p
        r = rnew
        j = j+1
        conv.append(np.linalg.norm(r)/np.linalg.norm(b))
    
        if b.all() == 0:
            err.append(np.linalg.norm(x))

    print(x, j, conv[j])


    plt.subplot(1, 2, 1)
    plt.plot(list(range(j+1)),conv, color = 'b')
    plt.xlabel('iteration')

    plt.subplot(1, 2, 2)
    plt.plot(list(range(j+1)),conv, color = 'b')
    plt.xlabel('iteration')
    plt.yscale('log')

    plt.show()


    return x 

A = np.array([[4,1], [1,3]])
print(A)

b = np.array([1,2])
print(b)

x_zero = np.array([2,1])
print(x_zero)

con_grad(A,x_zero,b)


