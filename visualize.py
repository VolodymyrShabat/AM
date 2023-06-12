import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import helpers


# def plot_AKT_cube(AKT,x_Quant,y_Quant,z_Quant):
#     x_coords = AKT[0]
#     y_coords = AKT[1]
#     z_coords = AKT[2]

#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')

#     # Scatter points
#     ax.scatter(x_coords, y_coords, z_coords, c='b', label='Points')
#     i = 0
#     j=0
#     while i<len(x_coords):
#         for k in range(y_Quant+1):
#             for j in range(i,i+(2*x_Quant)):
#                 if i>=len(x_coords)-1:
#                     break
#                 ax.plot([x_coords[i],x_coords[i+1]],[y_coords[i],y_coords[i+1]],[z_coords[i],z_coords[i+1]], c='gray')
#                 i+=1
#             i+=x_Quant+2
#         i+= (x_Quant+1)*(y_Quant)   
#         topPlus =(2*x_Quant+1)*(y_Quant+1)+(x_Quant+1)*y_Quant+(y_Quant+1)*(x_Quant+1)
#         topMiddlePlus = (2*x_Quant+1)*(y_Quant+1)+(x_Quant+1)*y_Quant
#         t=0
#         test=0
#         for yi in range(0,2*y_Quant+1):
#             print(t,"geing")
#             for xi in range(x_Quant+1):
#                 for _ in range(z_Quant):
#                     if t+topPlus>=len(x_coords):
#                         break
#                     ax.plot([x_coords[t+2*yi+yi],x_coords[t+topMiddlePlus-xi]],[y_coords[t],y_coords[t+topMiddlePlus-xi]],[z_coords[t],z_coords[t+topMiddlePlus-xi]], c='gray')
#                     t+=topMiddlePlus
#                     ax.plot([x_coords[t-xi],x_coords[t+(topPlus-topMiddlePlus)]],[y_coords[t-xi],y_coords[t+(topPlus-topMiddlePlus)]],[z_coords[t-xi],z_coords[t+(topPlus-topMiddlePlus)]], c='gray')
#                     t+=topPlus-topMiddlePlus
#                 t=(xi+1)*2
#             print(t,yi)
#         print("again")



#    # all nodes on top
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     ax.set_title('AKT Cube')
#     ax.legend()
#     plt.show()

class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __key(self):
        return self.x, self.y, self.z

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Point):
            return self.__key() == other.__key()

        return NotImplemented

def plot_cube(AKT, NT, elements):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    points = set()

    for element_index in range(len(NT[0])):
        visible_sides = get_visible_sides(element_index, elements)

        for side_index in range(1,7):
            if side_index not in visible_sides:
                continue

            side = np.empty([4, 9])

            for point_index in range(9):
                akt_index = NT[helpers.points_of_side[side_index][point_index]][element_index]
                pointX = AKT[0][akt_index]
                pointY = AKT[1][akt_index]
                pointZ = AKT[2][akt_index]
                side[0][point_index] = pointX
                side[1][point_index] = pointY
                side[2][point_index] = pointZ

            ax.plot(side[0], side[1], side[2], color='black')

            for point_index in range(8):
                    points.update([Point(
                        side[0][point_index],
                        side[1][point_index],
                        side[2][point_index],
                    )])

    for point in points:
        ax.scatter(point.x, point.y, point.z, color="c")

    ax.set_aspect('equal')
    plt.show()

def get_visible_sides(element_index, elements):
    visible_sides = set()

    for side_index, side in enumerate(elements[element_index]):
        if len(side) > 0:
            visible_sides.add(side_index)

    return visible_sides


def get_visible_sides(element_index, elements):
    x_pos, y_pos, z_pos = get_xyz_positions(element_index, elements)

    corners = get_corners(elements, x_pos, y_pos, z_pos)
    edges = get_edges(elements, x_pos, y_pos, z_pos)
    sides = get_sides(elements, x_pos, y_pos, z_pos)

    return corners | edges | sides


def get_xyz_positions(element_index, elements):
    count = 0

    for z in range(elements[2]):
        for y in range(elements[1]):
            for x in range(elements[0]):
                if count == element_index:
                    return x, y, z

                count += 1

    raise IndexError("Element index is out of bounds.")

def get_corners(elements, x_pos, y_pos, z_pos):
    result = set()

    if x_pos == 0 and y_pos == 0 and z_pos == 0:
        result.update([1, 4, 5])

    if x_pos == elements[0] - 1 and y_pos == 0 and z_pos == 0:
        result.update([1, 2, 5])

    if x_pos == 0 and y_pos == elements[1] - 1 and z_pos == 0:
        result.update([1, 3, 5])

    if x_pos == elements[0] - 1 and y_pos == elements[1] - 1 and z_pos == 0:
        result.update([2, 3, 5])

    if x_pos == 0 and y_pos == 0 and z_pos == elements[2] - 1:
        result.update([1, 4, 6])

    if x_pos == elements[0] - 1 and y_pos == 0 and z_pos == elements[2] - 1:
        result.update([1, 2, 6])

    if x_pos == 0 and y_pos == elements[1] - 1 and z_pos == elements[2] - 1:
        result.update([3, 4, 6])

    if x_pos == elements[0] - 1 and y_pos == elements[1] - 1 and z_pos == elements[2] - 1:
        result.update([4, 3, 6])

    return result

def get_edges(elements, x_pos, y_pos, z_pos):
    result = set()

    if y_pos == 0 and z_pos == 0:
        result.update([1, 5])

    if x_pos == 0 and z_pos == 0:
        result.update([4, 5])

    if x_pos == 0 and y_pos == 0:
        result.update([1, 4])

    if x_pos == elements[0] - 1 and z_pos == 0:
        result.update([2, 5])

    if y_pos == elements[1] - 1 and z_pos == 0:
        result.update([3, 5])

    if x_pos == elements[0] - 1 and y_pos == 0:
        result.update([1, 2])

    if x_pos == 0 and y_pos == elements[1] - 1:
        result.update([3, 4])

    if x_pos == elements[0] - 1 and y_pos == elements[1] - 1:
        result.update([2, 3])

    if x_pos == 0 and z_pos == elements[2] - 1:
        result.update([4, 6])

    if x_pos == elements[0] - 1 and z_pos == elements[2] - 1:
        result.update([2, 6])

    if y_pos == elements[1] - 1 and z_pos == elements[2] - 1:
        result.update([3, 6])

    if y_pos == 0 and z_pos == elements[2] - 1:
        result.update([1, 6])

    return result
      
def get_sides(elements, x_pos, y_pos, z_pos):
    if x_pos == 0:
        return {3}

    if y_pos == 0:
        return {0}

    if z_pos == 0:
        return {4}

    if x_pos == elements[0] - 1:
        return {1}

    if y_pos == elements[1] - 1:
        return {2}

    if z_pos == elements[2] - 1:
        return {5}

    return set()