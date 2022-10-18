import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from anim import Render, load_uv_data
import argparse
import numpy as np


if __name__ == "__main__":
    import sys
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str,
                        help="Folder containing u and v data", required=True)
    parser.add_argument('--cmap', type=str, default="jet", help="Colormap")
    parser.add_argument('--fps', type=int, default=15, help="FPS")
    parser.add_argument('--output', type=str, default="out.mp4", help="Output video")

    FLAGS, unparsed_args = parser.parse_known_args()
    if len(unparsed_args):
        print("Warning: unknow arguments {}".format(unparsed_args))

    metadata = dict(title='Movie Test', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps=FLAGS.fps, metadata=metadata)

    fig = plt.figure()
    fig.set_size_inches((8, 6))

    folder = FLAGS.folder
    #folder = "cpp-Release"
    udata, vdata = load_uv_data(folder)

    if len(udata) != len(vdata):
        print("Error with output data")
        sys.exit(1)

    cmap = plt.cm.get_cmap(FLAGS.cmap)
    fig = plt.figure()
    fig.set_size_inches((14, 6))
    ax1 = plt.subplot(1, 2, 1)
    imgU = plt.imshow(np.zeros_like(udata[0]), cmap=cmap)
    plt.colorbar(imgU).set_ticks([])
    ax1.set_title("U data", size='x-large')

    ax2 = plt.subplot(1, 2, 2)
    imgV = plt.imshow(np.zeros_like(vdata[0]), cmap=cmap)
    plt.colorbar(imgV).set_ticks([])
    ax2.set_title("V data", size='x-large')
    fig.tight_layout()

    render = Render(u=udata, v=vdata, imgU=imgU, imgV=imgV)
    with writer.saving(fig, FLAGS.output, 100):
        for i in range(len(udata)):
            render.update(i)
            writer.grab_frame()
            sys.stdout.write("frame -> {}\r".format(i+1))
            sys.stdout.flush()
    print()
