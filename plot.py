import numpy as np
import matplotlib
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams.update({'font.size': 12})

# Read the data file:
with open("out/data.txt", "r") as file:
    
    line = file.readline().split()
    n_grid = int(line[0])
    n_iter = int(line[1])
    
    line = file.readline().split()
    xmin = float(line[0])
    xmax = float(line[1])
    ymin = float(line[2])
    ymax = float(line[3])
    
    line = file.readline().split()
    data = [float(n) for n in line]


# Create a colormap:
Ncolors = 256
colors = [(0,'black'),(1/Ncolors,'midnightblue'),(0.2,'blue'),(0.35,'lightskyblue'),(0.5,'aliceblue'),(0.6,'lightyellow'),(0.75,'gold'),(0.85,'darkorange'),(0.95,'saddlebrown'),(1,'xkcd:chocolate brown')]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap', colors, N=Ncolors, gamma=1.0)

# Set the extent of the coordinate axes:
dx = (xmax-xmin)/(n_grid-1)
dy = (ymax-ymin)/(n_grid-1)
extent = [xmin-dx, xmax+dx, ymin-dy, ymax+dy]

# Plot the data:
X = (np.array(data)).reshape((n_grid,n_grid)).T
plt.imshow(X, origin='lower', extent=extent, cmap=cmap)
plt.xlabel('x')
plt.ylabel('y', rotation=0)

# set the color of the zoom window frame:
current_cmap = matplotlib.cm.get_cmap()
current_cmap.set_bad(color='white')
#plt.colorbar()

plt.show()
print(f"n_grid = {n_grid}, n_iter = {n_iter}")
