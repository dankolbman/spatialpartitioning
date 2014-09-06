import Base.show

# Type for particle representation
type Part
  id::Int
  pos::Array{Float64}
  vel::Array{Float64}
end

# Type for cell grids
type Cell
  id::Int
  parts::Array{Part}
  neighbors::Array{Cell}
end

# For printing the particles
function show(io::IO, c::Cell)
  print(io, "$(c.id)")
end

# For printing the particles
function show(io::IO, p::Part)
  print(io, "$(p.pos[1]) $(p.pos[2])\n")
end

# Runs a simulation with N particles for specified number of steps in a box
# of size S
function simulation(N=1000, steps=1000, S = 10, D=1)

  # Initialize system
  parts = initParts(N, S)
  cells = initCells(D)

  tic()
  @inbounds for step in 1:steps
    assignParts(parts, cells, S, D)
    collisionCheck(cells)
    #naiveCollision(parts)
    
    # Update positons
    @inbounds for p in parts
      npos = p.pos .+ p.vel
      # Keep in box
      if( npos[1] > 0.5 && npos[1] < S-0.5 && npos[2] > 0.5 && npos[2] < S-0.5)
        p.pos = npos
      end
      p.vel = [0.0, 0.0] #(2.*rand(2).-1.0)*0.4
    end
  end
  toc()
  println("DONE")
  writeParts(parts)
  run(`python posplot2D.py parts.dat`)
end

function naiveCollision(parts)
  @inbounds for part1 in 1:(size(parts,1)-1)
    @inbounds for part2 in 2:size(parts,1)
      p1 = parts[part1]
      p2 = parts[part2]
      dr = [(p1.pos[1] - p2.pos[1]) , (p1.pos[2] - p2.pos[2])]
      if(p1.id != p2.id)
        dr2 = dr[1]^2+dr[2]^2
        ang = atan2(dr[2], dr[1])
        rep = 0
        if( dr2 < (1.0))
          rep = -1*(1.0 - 1.0/sqrt(dr2))
        end
        p1.vel = p1.vel .+ [ rep*cos(ang), rep*sin(ang) ]
        p2.vel = p2.vel .- [ rep*cos(ang), rep*sin(ang) ]
      end
    end
  end
end

# Gets a cell at the given location after checking that it is in bounds
# Parameters:
#   cellGrid - The 2D grid of cells to get a cell from
#   row - The row coordinate
#   col - The col coordinate
#   D - The number of cells per dimension
# Returns:
#   A cell or None if request coord is out of bounds
function getNeighbor( cellGrid, row, col, D )
  if( row > 0 && col > 0 && row <= D && col <= D)
    return cellGrid[row, col]
  else
    return None
  end
end

# Initializes cells by instantiating them and linking them in a grid
# Params:
#   D - The number of grid tiles per dimension
# Returns:
#   A 2D Array of cells
function initCells(D)
  tCell = D*D
  cellGrid = Array(Cell, D, D)
  # Create cells
  for i in 1:tCell
    cell = Cell(i, Array(Cell,0), Array(Part,0))
    cellGrid[i] = cell
  end

  # Link cells
  for col in 1:D
    for row in 1:D
      n = Array(Cell,0)
      
      # Top
      if(getNeighbor(cellGrid, row+1, col, D) != None)
        push!(n, getNeighbor(cellGrid, row+1, col, D))
      end
      # Top Left
      if(getNeighbor(cellGrid, row+1, col-1, D) != None)
        push!(n, getNeighbor(cellGrid, row+1, col-1, D))
      end
      # Left
      if(getNeighbor(cellGrid, row, col-1, D) != None)
        push!(n, getNeighbor(cellGrid, row, col-1, D))
      end
      # Bottom Left
      if(getNeighbor(cellGrid, row-1, col-1, D) != None)
        push!(n, getNeighbor(cellGrid, row-1, col-1, D))
      end
      # Bottom
      if(getNeighbor(cellGrid, row-1, col, D) != None)
        push!(n, getNeighbor(cellGrid, row-1, col, D))
      end
      # Bottom Right
      if(getNeighbor(cellGrid, row-1, col+1, D) != None)
        push!(n, getNeighbor(cellGrid, row-1, col+1, D))
      end
      # Right
      if(getNeighbor(cellGrid, row, col+1, D) != None)
        push!(n, getNeighbor(cellGrid, row, col+1, D))
      end
      # Top Right
      if(getNeighbor(cellGrid, row+1, col+1, D) != None)
        push!(n, getNeighbor(cellGrid, row+1, col+1, D))
      end

      # Update neighbor list
      cellGrid[row,col].neighbors = n
    end
  end
  return cellGrid
end

# Assign particles to cells
function assignParts(parts, cellGrid, S, D)
  cSize = S/D
  # Clear all current particle lists for cells
  for c in cellGrid
    c.parts = Array(Part, 0)
  end
  for p in parts
    hash = div(p.pos[1], cSize) + div(p.pos[2], cSize)*D + 1
    push!(cellGrid[hash].parts, p)
  end
end

# Test for collisions
function collisionCheck(cellGrid)
  for c1 in cellGrid
    # Check particles within current cell
    for p1 in c1.parts
      for p2 in c1.parts
        collide(p1, p2)
      end
    end
    # Check collisions with neighboring cells
    for c2 in c1.neighbors
      collideCells(c1, c2)
    end
  end
end

# Collides all the particles within two cells together
# Params:
#   c1 - The first cell
#   c2 - The second cell
function collideCells(c1::Cell, c2::Cell)
  for p1 in c1.parts
    for p2 in c2.parts
      collide(p1, p2)
    end
  end
end

# Collides two particles together
# Params:
#   p1 - The first particle
#   p2 - The second particle
function collide(p1::Part, p2::Part)
  if(p1.id < p2.id)
    dr = [(p1.pos[1] - p2.pos[1]) , (p1.pos[2] - p2.pos[2])]
    dr2 = dr[1].*dr[1]+dr[2].*dr[2]
    ang = atan2(dr[2], dr[1])
    rep = 0
    if( dr2 < (1.0))
      rep = -1*(1.0 - 1.0/sqrt(dr2))
      p1.vel = p1.vel .+ [ rep*cos(ang), rep*sin(ang) ]
      p2.vel = p2.vel .- [ rep*cos(ang), rep*sin(ang) ]
    end
  end
end

# Creates particles in a box of size S
function initParts(N::Int64, S)
  # Holds particle objects
  parts = Array(Part, N)

  # Generate the particles
  parts = Array(Part, N)
  for p in 1:N
    pos = rand(2).*(S-2).+1.0
    parts[p] = Part(p, pos, [0.0,0.0])
  end

  return parts

end

# Write particles to file
function writeParts(parts)
  f = open("parts.dat", "w")
  for p in parts
    print(f, p)
  end
  close(f)
end

function main()
  # Number of particles
  N = 128000
  steps = 1000
  simulation(N, steps, 500, 500)
end

main()
