# Testing
if __name__ == "__main__":
    # Imports and setttings
    #import os
    #import sys
    #sys.path.insert(0, "../../"+os.getcwd())
    #print("add path.... ../../"+os.getcwd())
    
    import matplotlib.pyplot as plt
    from sklearn import datasets
    from core.LDA import *



    data = datasets.load_iris()
    X, y = data.data, data.target

    # Project the data onto the 2 primary linear discriminants
    lda = LDA(2)
    lda.fit(X, y)
    X_projected = lda.transform(X)

    print("Shape of X:", X.shape)
    print("Shape of transformed X:", X_projected.shape)

    x1, x2 = X_projected[:, 0], X_projected[:, 1]
    
    fig, ax = plt.subplots(2)
    im1 = ax[1].scatter(
        x1, x2, c=y, edgecolor="none", alpha=0.8, cmap=plt.cm.get_cmap("viridis", 3)
    )

    ax[1].set(xlabel="Linear Discriminant 1", ylabel="Linear Discriminant 2")
    fig.colorbar(im1,ax=ax[1], orientation='vertical')

    
    im2 = ax[0].scatter(
        X[:, 0], X[:, 1], c=y, edgecolor="none", alpha=0.8, cmap=plt.cm.get_cmap("viridis", 3)
    )
    ax[0].set(xlabel="Linear Discriminant 1", ylabel="Linear Discriminant 2")
    fig.colorbar(im1,ax=ax[0], orientation='vertical')
    
    plt.show()
