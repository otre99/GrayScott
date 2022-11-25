import numpy as np
from scipy.ndimage import laplace
import argparse
from matplotlib.colors import LightSource
import time

class GrayScottSover:
    def __init__(self, kparam: float, fparam: float, mu: float, mv: float, dim: int):
        self.kparam = kparam
        self.fparam = fparam
        self.mu = mu
        self.mv = mv
        self.u = np.zeros((dim, dim))
        self.v = np.zeros((dim, dim))
        self.uTemp = self.u.copy()
        self.vTemp = self.v.copy()

        self.dim = dim
        self.inter_count = 0

        r = self.dim // (self.dim // 16)
        r1 = self.dim // 2 - r
        r2 = self.dim // 2 + r

        self.u[...] = 1.0
        self.u[r1:r2, r1:r2] = 0.5
        self.v[...] = 0.0
        self.v[r1:r2, r1:r2] = 0.25

    def display_parameters(self):
        print("Parameters:")
        print("  Dims   -> {}x{}".format(self.dim, self.dim))
        print("  mu     -> {}".format(self.mu))
        print("  mv     -> {}".format(self.mv))
        print("  kparam -> {}".format(self.kparam))
        print("  fparam -> {}".format(self.fparam))

    def step(self):
        laplace(input=self.u, output=self.uTemp, mode='wrap')
        laplace(input=self.v, output=self.vTemp, mode='wrap')
        uv2 = self.u*self.v**2
        self.u += self.uTemp*self.mu - uv2 + self.fparam * (1-self.u)
        self.v += self.vTemp*self.mv + uv2 - (self.fparam+self.kparam)*self.v

        self.inter_count += 1


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


def create_cmd_parser():
    p = argparse.ArgumentParser(usage="GrayScott python implementation")
    p.add_argument("--fparam", default=0.055, type=float, help="f parameter")
    p.add_argument("--kparam", default=0.062, type=float, help="f parameter")
    p.add_argument("--mu", default=0.16, type=float, help="mu parameter")
    p.add_argument("--mv", default=0.08, type=float, help="mv parameter")
    p.add_argument("--dim", default=256, type=int, help="dim parameter")
    p.add_argument("--nsteps", default=10000, type=int, help="nsteps parameter")
    p.add_argument('--cmap', type=str, default="jet", help="Colormap")
    return p


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    parser = create_cmd_parser()
    FLAGS, UNKNOWN = parser.parse_known_args()

    if UNKNOWN:
        print("Warning: Unknown arguments {}".format(UNKNOWN))
    gs = GrayScottSover(fparam=FLAGS.fparam, kparam=FLAGS.kparam, mu=FLAGS.mu, mv=FLAGS.mv, dim=FLAGS.dim)
    gs.display_parameters()
    print("  nsteps -> {}".format(FLAGS.nsteps))

    t0 = time.time()
    for _ in range(FLAGS.nsteps):
        gs.step()
    t1 = time.time()
    print("Elapsed time {:.4f} seconds".format(t1-t0))

    fig = plt.figure()
    fig.set_size_inches((14, 6))
    ax = plt.subplot(1, 2, 1)
    show_data(ax=ax, data=gs.u, cmap=FLAGS.cmap, title="U data")
    ax = plt.subplot(1, 2, 2, sharex=ax, sharey=ax)
    show_data(ax=ax, data=gs.v, cmap=FLAGS.cmap, title="V data")
    plt.tight_layout()

    plt.show()
