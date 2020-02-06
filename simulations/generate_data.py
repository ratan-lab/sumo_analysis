from sklearn.datasets import make_blobs
from pandas import DataFrame
import logging
import numpy as np
import os
import argparse


def generate_clusters(nsamples: int, nfeatures: int, nclusters: int, cluster_std: int = 1.0, rstate: int = None, sd_dist: float = 1.0):
    """ Generate linearly separate clusters of data points
    Args:
        nsamples (int): number of points in all clusters
        nfeatures (int): number of features (at least 2)
        nclusters (int): number of clusters
        cluster_std (int): the standard deviation of the clusters
        rstate (int): optional, determines random number generation for dataset creation
        sd_dist (float): number of standard deviations between cluster centroids
    Returns:
        x (Numpy.ndarray): matrix of data point coordinates (features x samples)
        y (Numpy.ndarray): vector of data point labels (1 x samples)
    """
    assert all([x > 0 for x in [nsamples, nfeatures, nclusters]])

    # generate classification dataset
    centers = np.array([[i * cluster_std * sd_dist] * nfeatures for i in range(nclusters)])
    x, y = make_blobs(n_samples=nsamples, n_features=nfeatures, centers=centers, cluster_std=cluster_std,
                      random_state=rstate)

    # points grouped by their labels
    data = {"x" + str(i): x[:, i] for i in range(x.shape[1])}
    data["label"] = y
    df = DataFrame(data)
    grouped_x, grouped_y = [], []
    for key, group in df.groupby('label'):
        grouped_x.append(group.values[:, :-1])
        grouped_y.append(group.values[:, -1])

    return np.concatenate(tuple(grouped_x)), np.concatenate(tuple(grouped_y))


def generate_data(nsamples: int, nfeatures: int, nclusters: int, outdir: str, rstate: int = None, sd_dist: float = 1.0):
    # create output directory
    os.makedirs(outdir, exist_ok=True)

    # create data
    fmat, lab = generate_clusters(nsamples=nsamples, nfeatures=nfeatures, nclusters=nclusters, rstate=rstate, sd_dist=sd_dist)

    # save features
    fname = os.path.join(outdir, 'base_data.tsv')
    data = DataFrame(fmat.T, columns=['sample_{}'.format(x) for x in np.arange(nsamples)],
                     index=['feature_{}'.format(x) for x in np.arange(nfeatures)])
    data.to_csv(fname, sep="\t")
    logging.info("Data saved to: {}".format(fname))

    # save labels
    fname = os.path.join(outdir, 'base_labels.tsv')
    labels = np.array([['sample_{}'.format(x) for x in np.arange(nsamples)], [int(x) for x in lab]]).T
    labels_df = DataFrame(labels, columns=['sample', 'label'], index=None)
    labels_df.to_csv(fname, sep="\t", index=False)
    logging.info("Labels saved to: {}".format(fname))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate data for simulations")
    parser.add_argument('nsamples', type=int, help='number of samples')
    parser.add_argument('nfeatures', type=int, help='number of features')
    parser.add_argument('nclusters', type=int, help='number of clusters')
    parser.add_argument('outdir', type=str, help='output directory')
    parser.add_argument('-rstate', action='store', type=int, required=False, default=None,
                        help='random state')    
    parser.add_argument('-sd', action='store', type=float, required=False, default=1.0,
                        help='number of standard deviations between cluster centroids')
    args = parser.parse_args()

    # set up logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s"
    )
    logging.info("#Outdir: {}".format(args.outdir))

    generate_data(nsamples=args.nsamples, nfeatures=args.nfeatures, nclusters=args.nclusters, outdir=args.outdir,
                  rstate=args.rstate, sd_dist=args.sd)
