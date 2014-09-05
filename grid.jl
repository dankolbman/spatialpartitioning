
import Base.show

type Part
  id::Int
  pos::Array{Float64}
  vel::Array{Float64}
end

# For printing the particles
function show(io::IO, p::Part)
  print(io, "$(p.pos[1]) $(p.pos[2])\n")
end

# Runs a simulation with N particles for specified number of steps in a box
# of size S
function simulation(N=1000, steps=1000, S = 10, D=1)

  # Initialize system
  # Holds home cell ID and IDs of phantom cells
  #cellID = Array(Int, N, 4)
  # Holds object id and home and phantom cell types
  #objList = Array(Int, N)
  # Holds control bits to specify type of cell
  #ctrlBits = Array(Int, N, 4)
  objBins = zeros(Int, N,3)
  cells = zeros(Int, D^2)
  parts = init(N, S)
  #writeParts(parts)
  #run(`python posplot2D.py parts.dat`)

  tic()
  @inbounds for step in 1:steps
    # Update grid list
    objBins, cells = genGrid(parts, objBins, cells, S, D)
    #println(objBins)
    # We now have a sorted cell list ready for detection
    collisionCheck(parts, objBins, cells, D)
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

function hashPart(pos, cellSize, D)
  xhash =  div(pos[1],cellSize)
  yhash =  div(pos[2],cellSize)
  hash = yhash*D + xhash + 1
  #hash = (int(floor(pos[1]/cellSize)) << 3) | (int(floor(pos[2]/cellSize)) << 6)
  return hash
end

function hashType(pos, cellSize)

  return (div(pos[1],cellSize))%2+1+(div(pos[2],cellSize)%2)*2
end

# Generate grid cells
function genGrid(parts, objBins, cells, S, D)
  # Determine cell size
  cellSize = S/D
  # Total cells
  tCell = D^2
  for p in parts
    ID = p.id
    hash = hashPart(p.pos, cellSize, D)
    objBins[ID,1] = hash
    objBins[ID,2] = 999
    objBins[ID,3] = ID
  end
  newBins = Array(Int, size(objBins, 1), 3)
  perm = sortperm(objBins[:,1], alg=MergeSort)
  for i in 1:size(objBins, 1)
    newBins[i,:] = objBins[perm[i],:]
  end
  objBins = newBins

  # Index where the cells begin in the objCell
  cells = zeros(Int, tCell)
  bin = -1
  for i in 1:size(objBins, 1)
    if(objBins[i,1] != bin)
      bin = objBins[i,1]
      cells[bin] = i
      if(bin >= tCell) break end
    end
  end

  return objBins, cells

end

# Test for collisions
function collisionCheck(parts, objBin, cells, D)
  tCell = D*D
  N = size(parts,1)

  # Iterate each cell
  for i in 1:size(cells,1)
    # Check that there are paticles in the cell
    if( cells[i] > 0)
      # Where the cell begins
      j = cells[i]
      # The cell id of the current cell
      binj = objBin[j,1]
      # Iterate through the cell
      @inbounds while(j < N && objBin[j,1] == binj)
        binj = objBin[j,1]
        k = cells[i]
        @inbounds while( k <= N && objBin[k,1] == binj)
          collide(parts[objBin[j,3]], parts[objBin[k,3]])
          k+=1
        end
        # Iterate cell to the right
        right = i+1
        if(right <= tCell && cells[right] > 0)
          k = cells[right]
          bink = objBin[k]
          while( k < N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the top right
        tr = i + D + 1
        if(tr <= tCell && cells[tr] > 0)
          k = cells[tr]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the top
        top = i + D
        if(top <= tCell && cells[top] > 0)
          k = cells[top]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the top left
        tl = i + D - 1
        if(tl <= tCell && cells[tl] > 0)
          k = cells[tl]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the left
        left = i - 1
        if(left > 0 && cells[left] > 0)
          k = cells[left]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the bottom left
        bl = i - D - 1
        if(bl > 0 && cells[bl] > 0)
          k = cells[bl]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the bottom
        bottom = i - D
        if(bottom > 0 && cells[bottom] > 0)
          k = cells[bottom]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        # Iterate cell to the bottom right
        br = i - D + 1
        if(br > 0 && cells[br] > 0)
          k = cells[br]
          bink = objBin[k]
          @inbounds while(k <= N && objBin[k,1] == bink)
            collide(parts[objBin[j,3]], parts[objBin[k,3]])
            k+=1
          end
        end
        j+=1
      end
    end
  end

end

function collide(p1, p2)
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
function init(N::Int64, S)
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
  N = 32000
  steps = 100
  simulation(N, steps, 500, 500)

  #run(`python posplot2D.py parts.dat`)

end

main()
