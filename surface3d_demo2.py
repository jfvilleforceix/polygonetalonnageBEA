# -*- coding: utf-8 -*-
# Polygone d'étalonnage du BEA

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.tri as tri
from scipy.spatial import ConvexHull as CHu
from matplotlib import cm


# Write OBJ files from Zip object (X,Y,Z)
def writeOBJ(path_obj, zip, grid, geod):
    f = open(path_obj, "w")
    for xyz in zip:
        vertex = ' '.join(tuple(str(coord) for coord in xyz))
        vertex = "v " + vertex + "\n"
        f.write(vertex)
    if not geod:
        for simplex in grid:
            face = ' '.join(str(id+1) for id in simplex)
            face = "f " + face + "\n"
            f.write(face)
    doublons = []
    if geod:
        nblin = grid[0]
        nbcol = grid[1]
        print(grid)
        for i in range(nblin-1):
            for j in range(nbcol):
                no = i * nbcol + j
                ne = i * nbcol + j + 1
                so = (i + 1) * nbcol + j
                se = (i + 1) * nbcol + j + 1
                maille = [no, ne, so, se]
                if maille in doublons:
                    continue
                else:
                    doublons.append(maille)
                    face = ' '.join(str(id+1) for id in maille)
                    face = "f " + face + "\n"
                    f.write(face)
    f.close()


## Calibration polygon for 360 cameras


# Load icosahedral txt file
icosahedral_txt = open("ipack.3.110.txt","r")
icosahedral = [float(line[:-1]) for line in icosahedral_txt.readlines()]
icosahedral_txt.close()
# print(icosahedral)
x = []
y = []
z = []

for i in range(0,len(icosahedral),3):
    x.append(icosahedral[i])
    y.append(icosahedral[i+1])
    z.append(icosahedral[i+2])
# Export CSV and OBJ
rows = zip(x, y, z)
with open('icosaheral.csv', 'w', newline='') as icosaheral:
    csv_writer = csv.writer(icosaheral, delimiter=';')
    for row in rows:
        csv_writer.writerow(row)
rows = zip(x, y, z)
x = np.asarray(x, dtype=np.float).reshape(-1,1)
y = np.asarray(y, dtype=np.float).reshape(-1,1)
z = np.asarray(z, dtype=np.float).reshape(-1,1)
xyz = np.hstack((x,y,z))

# Keep only z-positive points
# xyzpos = np.asarray([row for row in xyz if row[2] >= 0])
xyzpos = np.asarray(xyz)

# Plot the icosahedral points cloud
# ax.plot_wireframe(x, y, z, color="r")
# ax.plot_trisurf(x, y, z, linewidth=0, antialiased=False)

# Triangulate icosahedron with convex hull
cvx = CHu(xyzpos)
tri = tri.Triangulation(xyzpos[:,0], xyzpos[:,1], triangles=cvx.simplices)

writeOBJ('icosaheral.obj', rows, cvx.simplices, False)

# Plot the triangulated surface
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# # Hide grid lines
# ax.grid(False)
# # plt.axis('off')
#
# ax.set_xlim(-1.01, 1.01)
# ax.set_ylim(-1.01, 1.01)
# ax.set_zlim(-0.01, 1.01)
# ax.scatter(xyzpos[:,0], xyzpos[:,1], xyzpos[:,2], color="r", s=20)
# ax.plot_trisurf(tri, xyzpos[:,2], color=(0,0,1,0), edgecolor='b', linewidth=0.5)
# ax.scatter(x, y, z, color='r')


## Calibration polygon for GoPro, standard cameras and zooms


# Function to create geodesic grid
def geodesicGrid(fov_hori, fov_verti, res_ang):
    demi_fov_verti_rad = fov_verti * np.pi / 360
    u = np.linspace(- fov_hori * np.pi / 360, fov_hori * np.pi / 360, int(np.ceil(fov_hori/res_ang)))
    v = np.linspace(np.pi/2 - fov_verti * np.pi / 180, np.pi/2, int(np.ceil(fov_verti/res_ang)))
    v += demi_fov_verti_rad
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    c, s = np.cos(-demi_fov_verti_rad), np.sin(-demi_fov_verti_rad)
    rot_y = np.array(((c, 0, s), (0, 1, 0), (-s, 0, c)))
    for i in range(int(np.ceil(fov_hori/res_ang))):
        for j in range(int(np.ceil(fov_verti/res_ang))):
            xyzij = [x[i][j], y[i][j], z[i][j]]
            # if rot_y.dot(xyzij)[2] >= 0:
            x[i][j] = rot_y.dot(xyzij)[0]
            y[i][j] = rot_y.dot(xyzij)[1]
            z[i][j] = rot_y.dot(xyzij)[2]
            # else:
            #     x[i][j] = 0
            #     y[i][j] = 0
            #     z[i][j] = 0
    return x, y, z


# GoPro
xGP, yGP, zGP = geodesicGrid(130, 80, 10)
rows = zip(np.concatenate(xGP, axis=0).tolist(), np.concatenate(yGP, axis=0).tolist(), np.concatenate(zGP, axis=0).tolist())
with open('goproGrid.csv', 'w', newline='') as goproGrid:
    csv_writer = csv.writer(goproGrid, delimiter=';')
    for row in rows:
        csv_writer.writerow(row)
rows = zip(np.concatenate(xGP, axis=0).tolist(), np.concatenate(yGP, axis=0).tolist(), np.concatenate(zGP, axis=0).tolist())
writeOBJ('goproGrid.obj', rows, xGP.shape, True)

# 50 mm
x50, y50, z50 = geodesicGrid(40, 27, 5)
rows = zip(np.concatenate(x50, axis=0).tolist(), np.concatenate(y50, axis=0).tolist(), np.concatenate(z50, axis=0).tolist())
with open('50mmGrid.csv', 'w', newline='') as a50mmGrid:
    csv_writer = csv.writer(a50mmGrid, delimiter=';')
    for row in rows:
        csv_writer.writerow(row)
rows = zip(np.concatenate(x50, axis=0).tolist(), np.concatenate(y50, axis=0).tolist(), np.concatenate(z50, axis=0).tolist())
writeOBJ('50mmGrid.obj', rows, x50.shape, True)

# 200 mm
x200, y200, z200 = geodesicGrid(11, 7, 2)
rows = zip(np.concatenate(x200, axis=0).tolist(), np.concatenate(y200, axis=0).tolist(), np.concatenate(z200, axis=0).tolist())
with open('200mmGrid.csv', 'w', newline='') as a200mmGrid:
    csv_writer = csv.writer(a200mmGrid, delimiter=';')
    for row in rows:
        csv_writer.writerow(row)
rows = zip(np.concatenate(x200, axis=0).tolist(), np.concatenate(y200, axis=0).tolist(), np.concatenate(z200, axis=0).tolist())
writeOBJ('200mmGrid.obj', rows, x200.shape, True)

# # Plot the geodesic grids
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# # Hide grid lines
# ax.grid(False)
# plt.axis('off')

ax.set_xlim(-1.01, 1.01)
ax.set_ylim(-1.01, 1.01)
ax.set_zlim(-0.01, 1.01)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
# ax.plot_surface(x, y, z, color='b')
# ax.scatter(xGP, yGP, zGP, color="r", s=10)
# ax.plot_wireframe(xGP, yGP, zGP, color="b", linewidth=0.5)
# ax.scatter(x50, y50, z50, color="r", s=10)
# ax.plot_wireframe(x50, y50, z50, color="b", linewidth=0.5)
ax.scatter(x200, y200, z200, color="r", s=10)
ax.plot_wireframe(x200, y200, z200, color="b", linewidth=0.5)
plt.show()

# triangulation a partir de points 3d

#
# rad = np.linalg.norm(xyzpos, axis=1)
# zen = np.arccos(xyzpos[:,2] / rad)
# azi = np.arctan2(xyzpos[:,1], xyzpos[:,0])

# tris = mtri.Triangulation(zen, azi)
#
# fig = plt.figure()
# # ax = fig.gca(projection='3d')
# ax  = fig.add_subplot(111, projection='3d')
# ax.plot_trisurf(xyzpos[:,0], xyzpos[:,1], xyzpos[:,2], triangles=tris.triangles)#, cmap=plt.cm.bone)
# plt.show()
