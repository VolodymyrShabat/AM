import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def plot_AKT_cube(AKT,x_Quant,y_Quant,z_Quant):
    x_coords = AKT[0]
    y_coords = AKT[1]
    z_coords = AKT[2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Scatter points
    ax.scatter(x_coords, y_coords, z_coords, c='b', label='Points')
    i = 0

    while i<len(x_coords):
        if i==(x_Quant+1)*(y_Quant+1)+(y_Quant)*(x_Quant+1)+(y_Quant+1)*x_Quant:
            print(i,"tjer")
            i+=5
            continue
        if i==4 and i!=0:
              print(i,"hello")
              i+=x_Quant+2
              print(i)
              continue
        if i<=len(x_coords)-3:
            print(i,"plot")
            ax.plot([x_coords[i],x_coords[i+1]],[y_coords[i],y_coords[i+1]],[z_coords[i],z_coords[i+1]], c='gray')
            i+=1
            continue

        
       
        
        i+=1

        
         
        
   

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('AKT Cube')
    ax.legend()
    plt.show()