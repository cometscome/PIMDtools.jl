module PIMDtools
    using LinearAlgebra
    export Trajectory,a0,get_relative_positions,calc_MSD
    const a0 = 0.5291772 #Bohr length [Angstrom]

    struct Atom
        x::Float64 #[a.u.]
        y::Float64
        z::Float64
        vx::Float64
        vy::Float64
        vz::Float64
        Fx::Float64
        Fy::Float64
        Fz::Float64
        E::Float64
        Atom(A) = new(A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10])
    end

    struct Snapshot
        atoms::Array{Atom,1}
        Snapshot(atoms) = new(atoms)
    end

    struct Trajectory
        data::Array{Snapshot,1}
        numtrj::Int64
        numatoms::Int64
        atomnames::Array{String,1}
        numeachatom::Array{Int64,1}
        periodic::Bool
        unitcell::Union{Nothing,Array{Array{Float64,1},1}}
        rotatematrix::Union{Nothing,Array{Float64,2}}
        rotatematrix_inv::Union{Nothing,Array{Float64,2}}
        itrj_start::Int64
    end

    Base.length(trj::Trajectory) = trj.numtrj

    #atominfo::Array{Tuple,1}
    #like: [("Pd",20),("Ru",16),("Al",92)] for Pd20Ru16Al92

    function Trajectory(filename,itrj_start,itrj_end,atominfo;
                            periodic=true,
                            unitcell=nothing
                        )

        numkinds = length(atominfo)
        print("System: ")
        atomnames = String[]
        numeachatom = Int64[]
        df = zeros(Float64,10)
        numtrj = itrj_end-itrj_start+1
        for (kind,num) in atominfo
            print("$(kind)$(num)")
            push!(atomnames,kind)
            push!(numeachatom,num)
        end
        println("\t")

        if periodic
            @assert unitcell != nothing "periocic system should have unitcell!"
            rotatematrix = zeros(Float64,3,3)
            for i=1:3
                rotatematrix[:,i] = unitcell[i]/norm(unitcell[i])
            end
            rotatematrix_inv = inv(rotatematrix)

            println("unitcell is ")
            for i=1:3
                println(unitcell[i]," [a.u.]")
            end
            println("rotatematrix is ",rotatematrix)
        else
            unitcell = nothing
            rotatematrix = nothing
            rotatematrix_inv = nothing
        end

        

        numatoms = sum(numeachatom)
        println("number of total atoms: ",numatoms)
        println("number of total data: ",numtrj)

        data = Array{Snapshot,1}(undef,numtrj)
        atoms = Array{Atom,1}(undef,numatoms)


        fp = open(filename,"r")
        for itrj=1:itrj_start-1
            for iatom=1:numatoms
                readline(fp)
            end
        end

        count = 0
        for itrj=itrj_start:itrj_end
            count += 1
            for iatom=1:numatoms
                u = split(readline(fp))
                di = parse(Int64,u[1]) 
                for i=2:11
                    df[i-1] = parse(Float64,u[i]) 
                end
                atoms[iatom] = Atom(df[1:10])
            end
            data[count] = deepcopy(Snapshot(atoms))
        end

        return Trajectory(
            data,
            numtrj,
            numatoms,
            atomnames,
            numeachatom,
            periodic,
            unitcell,
            rotatematrix,
            rotatematrix_inv,
            itrj_start
        )
        #=
            data::Array{Snapshot,1}
            numtrj::Int64
            numatoms::Int64
            atomnames::Array{String,1}
            numeachatom::Array{Int64,1}
            periodic::Bool
            unitcell::Union{Nothing,Array{Array{Float64,1},1}}
            rotatematrix::Union{Nothing,Array{Float64,2}}
            rotatematrix_inv::Union{Nothing,Array{Float64,2}}
            itrj_start::Int64
        =#
    end

    function get_relative_positions(trj::Trajectory)
        positions = get_positions(trj)
        centerofmass = calc_centerofmass(trj,positions)
        for itrj=1:length(trj)
            for iatom=1:trj.numatoms
                positions[:,iatom,itrj] -= centerofmass[:,itrj]
            end
        end
        return positions
    end

    function calc_Galilei_frame!(trj::Trajectory,positions)
        centerofmass = calc_centerofmass(trj,positions)
        for itrj=1:length(trj)
            for iatom=1:trj.numatoms
                positions[:,iatom,itrj] -= centerofmass[:,itrj]
            end
        end
    end

    function calc_centerofmass(trj::Trajectory,positions )
        centerofmass = zeros(Float64,3,length(trj))
        for itrj=1:length(trj)
            for iatom=1:trj.numatoms
                for i=1:3
                    centerofmass[i,itrj] += positions[i,iatom,itrj]/trj.numatoms
                end
            end
        end
        return centerofmass
    end

    function get_positions(trj::Trajectory)
        positions = zeros(Float64,3,trj.numatoms,length(trj))
        dr = zeros(Float64,3)
        tempvec = zeros(Float64,3)
        dur = zeros(Float64,3)
        unorm = zeros(Float64,3)
        for i=1:3
            unorm[i] = norm(trj.unitcell[i])
        end

        for itrj=1:length(trj)
            for iatom=1:trj.numatoms
                #println(itrj,"\t",iatom)
                positions[1,iatom,itrj] = trj[itrj][iatom].x
                positions[2,iatom,itrj] = trj[itrj][iatom].y
                positions[3,iatom,itrj] = trj[itrj][iatom].z
                if itrj > 1
                    dx = positions[1,iatom,itrj]-positions[1,iatom,itrj-1]
                    dy = positions[2,iatom,itrj]-positions[2,iatom,itrj-1]
                    dz = positions[3,iatom,itrj]-positions[3,iatom,itrj-1]

                    dr[1] = dx
                    dr[2] = dy
                    dr[3] = dz
                    mul!(dur,trj.rotatematrix_inv,dr)

                    for i=1:3
                        tempvec .= 0
                        tempvec[i] = ifelse(dur[i] > unorm[i]/3,-unorm[i],0)
                        mul!(dr,trj.rotatematrix,tempvec)
                        positions[:,iatom,itrj] .+= dr    

                        tempvec .= 0
                        tempvec[i] = ifelse(dur[i] < -unorm[i]/3,unorm[i],0)
                        mul!(dr,trj.rotatematrix,tempvec)
                        positions[:,iatom,itrj] .+= dr    
                    end

                end
            end
        end
        return positions
    end

    function Base.getindex(x::Snapshot,i)
        return x.atoms[i]
    end

    function Base.getindex(x::Trajectory,i)
        return x.data[i]
    end

    function Base.getindex(trj::Trajectory,iatom,itrj)
        #println(trj[i][j].x)
        return trj[itrj][iatom]
    end

    function calc_MSD(trj::Trajectory,dirname,period,start_itrj,maxsteps)
        positions = get_relative_positions(trj)
        calc_MSD(trj,dirname,positions,period,start_itrj,maxsteps)
    end

    function calc_MSD(trj::Trajectory,dirname,positions,period,start_itrj,maxsteps)
        msd = zeros(maxsteps)
        phi = zeros(trj.numatoms,maxsteps)
        ri = zeros(3)
        drj = zeros(3)
        dphi = 0
        count = 0
        istart = start_itrj - trj.itrj_start + 1
        @assert istart > 0 "no data available!"
        for itrj=istart:period:length(trj)
            count += 1
            for iatom=1:trj.numatoms
                for i=1:3
                    ri[i] = positions[i,iatom,itrj]*a0
                end
                for it = 1:maxsteps
                    jtrj = itrj+it-1
                    if itrj+maxsteps-1 > length(trj)
                        break
                    end
                    dphi = 0.0
                    for i=1:3
                        dphi += (positions[i,iatom,jtrj]*a0 - ri[i])^2
                        #drj[i] = positions[i,iatom,jtrj]*a0 - ri[i]
                    end
                    msd[it] += dphi
                    phi[iatom,it] += dphi
                end

            end
            fp = open(dirname*"/atomic_MSD_from$(start_itrj)-th.txt","w")
            fp2 = open(dirname*"/MSD_from$(start_itrj)-th.txt","w")


            for it = 1:maxsteps
                println(fp2,it,"\t",msd[it]/(count*trj.numatoms))
                print(fp,it,"\t")
                for iatom=1:trj.numatoms
                    print(fp,phi[iatom,it]/count,"\t")
                end
                println(fp,"\t")
            end
            close(fp)
    
            close(fp2)
            println("itrj = ",itrj)
        end
        return msd,phi
    end


# Write your package code here.

end
