"""
  Plot positions of a multispecies system in 2d

  Dan Kolbman 2014
"""
import matplotlib.pyplot as plt
import sys

def readParts(filen):
  xpos = []
  ypos = []
  with open(filen,'r') as f:
    for line in f:
      l = line.split()
      xpos.append(l[0])
      ypos.append(l[1])
  return xpos, ypos

def plotSys():
  """
  Plots position data
  """
  fig = plt.gcf()
  ax = fig.gca()
  xpos, ypos = readParts('parts.dat')
  circleScatter(xpos, ypos, ax, radius=0.5)
  #s = [conf['diameter']**2/4*3.1415 for i in range(len(xpos))]
  #plt.scatter(xpos, ypos,color=colors[(i)%3])

  plotBounds( plt.gcf().gca())
  plt.savefig('fpos.png', figsize=(1,1), dpi=100)
  plt.grid()


def plotBounds( axes):
  """ plotBounds : Dict Axes -> True
  Makes a rectangle or circle depending on the geometry of the boundary
  """
  plt.xlim((-10, 510))
  plt.ylim((-10, 510))
  shape = plt.Rectangle((0, 0), 510, 510)
  shape.fill = False
  axes.add_artist(shape)
  return True
  
def circleScatter(xpos, ypos, axes, **kwargs):
  """ circleScatter : float[] float[] -> True
  Creates a scatter plot of circles
  """
  #for x,y in zip(xpos, ypos):
  for i in range(len(xpos)):
    #circle = plt.Circle((xpos[i],ypos[i]), color=str(i/len(xpos)), **kwargs)
    circle = plt.Circle((xpos[i],ypos[i]), color='b', **kwargs)
    axes.add_patch(circle)

  return True


"""
  If called directly, only show the position plot
"""
if __name__ == '__main__':
  plt.gcf().add_subplot(111, aspect='equal')
  plotSys()
  #plt.savefig(sys.argv[2])
  plt.show()
  

