import numpy as np
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scipy.interpolate import griddata
from matplotlib import pyplot as plt
from matplotlib import cm

def plot_3d_and_contour(filename, column, output_filename=None):
    data = np.loadtxt(filename)
    psi = data[:, 0]
    khi = data[:, 1]
    vals = data[:, int(column)]

    cols = vals.size
    # print("Number of points = %d" %cols)

    fig = plt.figure()
    # 3D plot of data
    ax = fig.add_subplot(211, projection='3d')
    surf = ax.plot_trisurf(psi, khi, vals, cmap=cm.jet, linewidth=0)
    fig.colorbar(surf)
    ax.set_xlim(psi.min(), psi.max())
    ax.set_ylim(khi.min(), khi.max())

    # Contour plot of data
    fig.add_subplot(212)
    # define grid.
    xi = np.linspace(psi.min(), psi.max(), 181)
    yi = np.linspace(khi.min(), khi.max(), 181)
    # grid the data.
    zi = griddata((psi, khi), vals, (xi[None, :], yi[:, None]), method='cubic')
    # print("zi = ", zi)
    # contour the gridded data
    CS = plt.contour(xi, yi, zi, 20, linewidths=0.5, colors='k')
    CS = plt.contourf(xi, yi, zi, 20, cmap=cm.jet)
    plt.colorbar()  # draw colorbar
    # plt.title('contour')
    # plt.show()
    if output_filename is None:
        fig.savefig(filename + '_col(' + str(column) + ').png')
    else:
        fig.savefig(output_filename)

    plt.close(fig)


def plot_3d_with_contour(filename, column, rownum=61, colnum=31, output_filename=None):
    data = np.loadtxt(filename)
    psi = data[:, 0]
    khi = data[:, 1]
    vals = data[:, 2]

    num_points = vals.size
    print("Number of points = %d\n" % num_points)

    assert (rownum * colnum == num_points), \
        "Number of rows and columns did not match with input data length!"

    X = psi.reshape(rownum, colnum)
    Y = khi.reshape(rownum, colnum)
    Z = vals.reshape(rownum, colnum)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, Z, rstride=4, cstride=2, alpha=0.3, linewidth=0.5, antialiased=True)
    # fig.colorbar(surf,shrink=0.75,aspect=5)
    cset = ax.contourf(X, Y, Z, 30, zdir='z', offset=Z.min(), cmap=cm.jet)
    # cset = ax.contourf(X, Y, Z, zdir='x', cmap=cm.coolwarm)
    # cset = ax.contourf(X, Y, Z, 30, zdir='y', offset = 45, cmap=cm.jet)

    ax.set_xlabel('$\psi$')
    # ax.set_xlim(-45, 45)
    ax.set_ylabel('$\chi$')
    # ax.set_ylim(-40, 40)
    # ax.set_zlabel('FD',labelpad=20, verticalalignment = 'top')
    # ax.set_zlim(2.5, 3)

    # plt.show()
    if output_filename is None:
        fig.savefig(filename + '_col(' + str(column) + ').png', dpi=300)
    else:
        fig.savefig(output_filename)

    plt.close(fig)


if __name__ == "__main__":
    import sys
    plot_3d_and_contour(sys.argv[1], int(sys.argv[2]))
    #plot_3d_with_contour(sys.argv[1], int(sys.argv[2]))
