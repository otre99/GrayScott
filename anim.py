from email.policy import default
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os 
from plot import load_data
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable

def load_uv_data(folder):
    ufiles=[]
    vfiles=[]
    flist = os.listdir(folder)
    for fname in flist:
        if fname.startswith("UITER"):
            ufiles.append(os.path.join(folder, fname))
        if fname.startswith("VITER"):
            vfiles.append(os.path.join(folder, fname))
    ufiles.sort()
    vfiles.sort()            
    return [load_data(s) for s in ufiles], [load_data(s) for s in vfiles]

class Render:
    def __init__(self, u, v, imgU, imgV) -> None:
        self.u = u
        self.v = v 
        self.imgU, self.imgV = imgU, imgV

    def update(self, i):
        arr = self.u[i] 
        self.imgU.set_data(arr)
        self.imgU.set_clim(np.min(arr), np.max(arr))

        arr = self.v[i]  
        self.imgV.set_data(arr)
        self.imgV.set_clim(np.min(arr), np.max(arr))
        return self.imgU, self.imgV

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, help="Folder containing u and v data", required=True)
    parser.add_argument('--cmap', type=str, default="jet", help="Colormap")
    parser.add_argument('--delay', type=int, default=1, help="Delay")

    FLAGS, unparsed_args = parser.parse_known_args()
    if len(unparsed_args):
        print("Warning: unknow arguments {}".format(unparsed_args))

    import sys 
    folder = FLAGS.folder
    udata, vdata = load_uv_data(folder)

    if len(udata)!=len(vdata):
        print("Error with output data")
        sys.exit(1)
    
    cmap=plt.cm.get_cmap(FLAGS.cmap)
    fig = plt.figure()
    ax1 = plt.subplot(1,2,1)
    imgU = plt.imshow(np.zeros_like(udata[0]),cmap=cmap)
    plt.colorbar(imgU).set_ticks([])
    ax1.set_title("U data", size='x-large')


    ax2 = plt.subplot(1,2,2)
    imgV = plt.imshow(np.zeros_like(vdata[0]),cmap=cmap)
    plt.colorbar(imgV).set_ticks([])
    ax2.set_title("V data", size='x-large')


    render = Render(u=udata, v=vdata, imgU=imgU, imgV=imgV)
    ani = FuncAnimation(plt.gcf(), render.update, frames=np.arange(len(udata)), blit=True, interval=FLAGS.delay)
    plt.show()
