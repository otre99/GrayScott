import argparse
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from matplotlib.colors import LightSource, Normalize


def show_data(ax, data, cmap, title=""):
    cmap = plt.cm.get_cmap(cmap)
    ls = LightSource(315, 45)
    rgb = ls.shade(data, cmap)
    ax.imshow(rgb, interpolation='bilinear')
    # Use a proxy artist for the colorbar...
    im = ax.imshow(data, cmap=cmap)
    im.remove()

    plt.colorbar(im, ax=ax)
    ax.set_title(title, size='x-large')


def load_data(fname):
    data = np.fromfile(file=fname, dtype=np.float64)
    n = int(sqrt(data.size))
    return data.reshape(n, n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--u', type=str, help="U data")
    parser.add_argument('--v', type=str, help="V data")
    parser.add_argument('--cmap', type=str, default="jet", help="Colormap")
    FLAGS, unparsed_args = parser.parse_known_args()
    if len(unparsed_args):
        print("Warning: unknow arguments {}".format(unparsed_args))

    u = load_data(fname=FLAGS.u) if FLAGS.u is not None else None
    v = load_data(fname=FLAGS.v) if FLAGS.v is not None else None

    fig = plt.figure()
    fig.set_size_inches((14, 6))
    if u is not None and v is not None:
        ax = plt.subplot(1, 2, 1)
        show_data(ax=ax, data=u, cmap=FLAGS.cmap, title="U data")
        ax = plt.subplot(1, 2, 2, sharex=ax, sharey=ax)
        show_data(ax=ax, data=v, cmap=FLAGS.cmap, title="V data")
    elif FLAGS.u is not None:
        show_data(ax=plt.gca(), data=u, cmap=FLAGS.cmap, title="U data")
    elif FLAGS.v is not None:
        show_data(ax=plt.gca(), data=v, cmap=FLAGS.cmap, title="V data")
    fig.tight_layout()
    plt.show()
