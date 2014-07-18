
import Base.show

type Part
  pos::Array{Float64}
  vel::Array{Float64}
end

# For printing the particles
function show(io::IO, p::Part)
  print(io, "$(p.pos[1]) $(p.pos[2])\n")
end

# Runs a simulation with N particles for specified number of steps in a box
# of size S
function simulation(N=1000, steps=1000, S = 10)

  # Initialize system
  parts = init(N, S)

  for step in 1:steps
    # Update grid list
    grid = getGrid(parts,S)
    # Update positons
    for p in parts
      p.vel = 2.*rand(2).-1.0
      npos = p.pos + p.vel
      # Keep in box
      if( npos[1] > 0.5 && npos[1] < S-0.5 && npos[2] > 0.5 && npos[2] < S-0.5)
        p.pos = npos
      end
    end
  end
  writeParts(parts)
end

function hashPart(p, invCellSize, D)
  xhash =  int(p.pos[1]*invCellSize)
  yhash =  int(p.pos[2]*invCellSize)
  println(yhash*D + xhash + 1)
  return yhash*D + xhash + 1
end

# Generate grid cells
function getGrid(parts, S)
  # Partition axis into D cells
  D = 5
  # Determine inve cell size
  invCellSize = D/S
  # Total cells
  totCells = D^2
  grid = Array(Any,totCells)
  println(size(grid,1))
  for cell in 1:size(grid,1)
    grid[cell] = Array{Part}
  end
  println(grid)
  for p in parts
     push!(grid[hashPart(p, invCellSize,D)], p)
  end
  println(grid)

end

# Creates particles in a box of size S
function init(N::Int64, S)

  parts = Array(Part, N)
  for p in 1:N
    pos = rand(2).*(S-2).+1.0
    parts[p] = Part(pos, [0,0])
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
  N = 10
  steps = 1
  simulation(N, steps, 100)

  run(`python posplot2D.py`)

end

main()
