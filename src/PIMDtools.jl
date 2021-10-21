module PIMDtools
    using LinearAlgebra
    #using HDF5
    export Trajectory,a0,get_relative_positions,calc_MSD,toBinary_r,
            calc_MSD_fromfile
    const a0 = 0.5291772 #Bohr length [Angstrom]

    const nbytes_Int64= sizeof(Int64)
    const nbytes_Float64= sizeof(Float64)

    struct Atom_vFE
        vx::Float64
        vy::Float64
        vz::Float64
        Fx::Float64
        Fy::Float64
        Fz::Float64
        E::Float64
    end

    struct Atom
        x::Float64 #[a.u.]
        y::Float64
        z::Float64
        vFE::Union{Nothing,Atom_vFE}

        Atom(A) = Atom(A,ronly=true)
        
        function Atom(A;ronly=true)
            if ronly
                return new(A[1],A[2],A[3],nothing)
            else
                return new(A[1],A[2],A[3],Atom_vFE(A[4],A[5],A[6],A[7],A[8],A[9],A[10]))
            end
        end
        #Atom(A) = new(A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10])
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

    Base.length(snapshot::Snapshot) = length(snapshot.atoms)

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

    #atominfo::Array{Tuple,1}
    #like: [("Pd",20),("Ru",16),("Al",92)] for Pd20Ru16Al92

    function Trajectory(filename,itrj_start,itrj_end,atominfo;
                            periodic=true,
                            unitcell=nothing,
                            ronly=true
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
                atoms[iatom] = Atom(df[1:10],ronly=ronly)
            end
            data[count] = deepcopy(Snapshot(atoms))
        end
        close(fp)

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

    """
        Binary:
        num_of_eachdata::Int64
        for itrj=start_itrj:start_itrj+num_of_eachdata-1
            itrj::Int64
            numatoms::Int64
            for iatom=1:numatoms
                x::Float64
                y::Float64
                z::Float64
            end
        end
    """
    function toBinary_r(filename,dirname,itrj_start,itrj_end,num_of_structures,atominfo)
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

        numatoms = sum(numeachatom)
        println("number of total atoms: ",numatoms)
        println("number of total data: ",numtrj)

        fp = open(filename,"r")
        for itrj=1:itrj_start-1
            for iatom=1:numatoms
                readline(fp)
            end
        end

        count = 0
        setcount = itrj_start
        countstruct = 0
        setstart = setcount
        setend  =setcount +num_of_structures - 1

        
        dataname = dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat"
        #ifp1 = h5open(dataname, "w")
        ifp1 = open(dataname, "w")
        num_of_eachdata = setend -setstart + 1
        write(ifp1,num_of_eachdata)
        #ifp1 = open(dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat","w")
        for itrj=itrj_start:itrj_end
            count += 1
            countstruct += 1
            write(ifp1,itrj)
            write(ifp1,numatoms)
            for iatom=1:numatoms
                u = split(readline(fp))
                di = parse(Int64,u[1]) 
                for i=2:11
                    df[i-1] = parse(Float64,u[i]) 
                end
                #atoms[iatom] = Atom(df[1:10],ronly=ronly)
                write(ifp1,df[1])
                write(ifp1,df[2])
                write(ifp1,df[3])
                #write(ifp1,"$(lpad(itrj,8,"0"))/$(lpad(iatom,4,"0"))/r",Float64[df[1],df[2],df[3]])
            end
            setcount += 1
            if countstruct == num_of_structures
                println("$dataname is done. $num_of_eachdata strctures")
                close(ifp1)                
                setstart = setcount
                
                if setstart <= itrj_end
                    setend  =setcount +num_of_structures - 1
                    if setend > itrj_end
                        setend =itrj_end
                    end
                    num_of_eachdata = setend -setstart + 1
                    
                    dataname = dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat"
                    #ifp1 = h5open(dataname, "w")
                    ifp1 = open(dataname, "w")
                    write(ifp1,num_of_eachdata)
                    countstruct = 0
                end
            end
        end
        println("$dataname is done. $num_of_eachdata strctures")
        close(ifp1) 
        close(fp)

    end

    function toBinary_r(filename,dirname,itrj_start,itrj_end,num_of_structures,atominfo,
        unitcell
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


        unorm = zeros(Float64,3)
        dr = zeros(Float64,3)
        tempvec = zeros(Float64,3)
        dur = zeros(Float64,3)
        for i=1:3
            unorm[i] = norm(unitcell[i])
        end



        numatoms = sum(numeachatom)
        println("number of total atoms: ",numatoms)
        println("number of total data: ",numtrj)

        fp = open(filename,"r")
        for itrj=1:itrj_start-1
            for iatom=1:numatoms
                readline(fp)
            end
        end

        count = 0
        setcount = itrj_start
        countstruct = 0
        setstart = setcount
        setend  =setcount +num_of_structures - 1

        df_struct = zeros(Float64,3,numatoms)
        df_struct_old = zeros(Float64,3,numatoms)
        dataname = dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat"
        #ifp1 = h5open(dataname, "w")
        ifp1 = open(dataname, "w")
        num_of_eachdata = setend -setstart + 1
        write(ifp1,num_of_eachdata)
        #ifp1 = open(dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat","w")
        for itrj=itrj_start:itrj_end
            count += 1
            countstruct += 1
            write(ifp1,itrj)
            write(ifp1,numatoms)
            for iatom=1:numatoms
                u = split(readline(fp))
                di = parse(Int64,u[1]) 
                for i=2:11
                    df[i-1] = parse(Float64,u[i]) 
                end
                for i=1:3
                    df_struct[i,iatom] = df[i]
                end

                if itrj > 1
                    dx = df_struct[1,iatom]-df_struct_old[1,iatom]
                    dy = df_struct[2,iatom]-df_struct_old[2,iatom]
                    dz = df_struct[3,iatom]-df_struct_old[3,iatom]

                    dr[1] = dx
                    dr[2] = dy
                    dr[3] = dz
                    mul!(dur,rotatematrix_inv,dr)

                    for i=1:3
                        tempvec .= 0
                        tempvec[i] = ifelse(dur[i] > unorm[i]/3,-unorm[i],0)
                        mul!(dr,rotatematrix,tempvec)
                        df_struct[:,iatom] .+= dr    

                        tempvec .= 0
                        tempvec[i] = ifelse(dur[i] < -unorm[i]/3,unorm[i],0)
                        mul!(dr,rotatematrix,tempvec)
                        df_struct[:,iatom] .+= dr    
                    end
                end
            end

            df_struct_old[:,:] = df_struct

            for iatom=1:numatoms
                #atoms[iatom] = Atom(df[1:10],ronly=ronly)
                for i=1:3
                    write(ifp1,df_struct[i,iatom])
                    #df_struct[i,iatom] = df[i]
                end
                #write(ifp1,"$(lpad(itrj,8,"0"))/$(lpad(iatom,4,"0"))/r",Float64[df[1],df[2],df[3]])
            end
            setcount += 1
            if countstruct == num_of_structures
                println("$dataname is done. $num_of_eachdata strctures")
                close(ifp1)                
                setstart = setcount
                
                if setstart <= itrj_end
                    setend  =setcount +num_of_structures - 1
                    if setend > itrj_end
                        setend =itrj_end
                    end
                    num_of_eachdata = setend -setstart + 1
                    
                    dataname = dirname*"/$(lpad(setstart,8,"0"))_$(lpad(setend,8,"0"))_trjset_r.dat"
                    #ifp1 = h5open(dataname, "w")
                    ifp1 = open(dataname, "w")
                    write(ifp1,num_of_eachdata)
                    countstruct = 0
                end
            end
        end
        println("$dataname is done. $num_of_eachdata strctures")
        close(ifp1) 
        close(fp)

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

    function get_relative_positions(positions_old,snapshot::Snapshot,
        unitcell,rotatematrix,rotatematrix_inv)
        positions,numatoms = get_positions(positions_old,snapshot,
                        unitcell,rotatematrix,rotatematrix_inv)
        centerofmass = calc_centerofmass(positions,numatoms)
        for iatom=1:numatoms
            for i=1:3
                positions[i,iatom] -= centerofmass[i]
            end
        end
        return positions
    end

    function get_relative_positions(snapshot::Snapshot)
        numatoms = length(snapshot)
        positions = zeros(Float64,3,numatoms)

        for iatom=1:numatoms
            positions[1,iatom] = snapshot[iatom].x
            positions[2,iatom] = snapshot[iatom].y
            positions[3,iatom] = snapshot[iatom].z
        end

        centerofmass = calc_centerofmass(positions,numatoms)
        for iatom=1:numatoms
            for i=1:3
                positions[i,iatom] -= centerofmass[i]
            end
        end
        return positions
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

    function calc_centerofmass(positions,numatoms)
        centerofmass = zeros(Float64,3)
        for iatom=1:numatoms
            for i=1:3
                centerofmass[i] += positions[i,iatom]/numatoms
            end
        end
        return centerofmass
    end

    #=
    function get_positions(positions_old,snapshot::Snapshot,
        unitcell,rotatematrix,rotatematrix_inv)
        
        unorm = zeros(Float64,3)
        dr = zeros(Float64,3)
        tempvec = zeros(Float64,3)
        dur = zeros(Float64,3)
        for i=1:3
            unorm[i] = norm(unitcell[i])
        end

        numatoms = length(snapshot)
        positions = zeros(Float64,3,numatoms)

        for iatom=1:numatoms
            positions[1,iatom] = snapshot[iatom].x
            positions[2,iatom] = snapshot[iatom].y
            positions[3,iatom] = snapshot[iatom].z
            dx = snapshot[iatom].x - positions_old[1,iatom]
            dy = snapshot[iatom].y - positions_old[2,iatom]
            dz = snapshot[iatom].z - positions_old[3,iatom]

            dr[1] = dx
            dr[2] = dy
            dr[3] = dz
            mul!(dur,rotatematrix_inv,dr)

            for i=1:3
                tempvec .= 0
                tempvec[i] = ifelse(dur[i] > unorm[i]/3,-unorm[i],0)
                mul!(dr,rotatematrix,tempvec)
                positions[:,iatom] .+= dr    

                tempvec .= 0
                tempvec[i] = ifelse(dur[i] < -unorm[i]/3,unorm[i],0)
                mul!(dr,rotatematrix,tempvec)
                positions[:,iatom] .+= dr    
            end
        end
        return positions,numatoms
    end
    =#

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



    struct Trjfile
        inputdirname::String
        filelist::Array{String,1}
        numdata_list::Array{Int64,1}
        index_list::Array{Int64,1}
        numatoms_list::Array{Int64,1}

        function Trjfile(inputdirname)
            filelist = filter(x -> x[end-3:end] == ".dat",readdir(inputdirname))
            println("filelist: ",filelist)
            numdata_list = Int64[]
            index_list = Int64[]
            numatoms_list = Int64[]
            count = 0
            for file in filelist
                println(file)
                fp = open(inputdirname*"/"*file,"r")
                #println(u)
                push!(numdata_list,read(fp,Int64))
                push!(index_list,read(fp,Int64))
                push!(numatoms_list,read(fp,Int64))
                seek(fp,0)
                close(fp)
            end

            return new(inputdirname,filelist,numdata_list,index_list,numatoms_list)
        end
    end

    function get_lastindex(trjfile::Trjfile)
        return trjfile.index_list[end] + trjfile.numdata_list[end] -1
    end

    function get_structure(trjfile::Trjfile,itrj)
        #fileindex = searchsortedfirst(trjfile.index_list,itrj)
        #println(filter(x -> x <= itrj,trjfile.index_list)[end])
        istart = filter(x -> x <= itrj,trjfile.index_list)[end]
        fileindex = searchsortedfirst(trjfile.index_list,istart)
        fp = open(trjfile.inputdirname*"/"*trjfile.filelist[fileindex],"r")

        seekposition = 0
        count = 0
        seekposition += nbytes_Int64
        println(istart)
        for i=istart:itrj-1
            count += 1
            seekposition += 2*nbytes_Int64
            for iatom=1:trjfile.numatoms_list[count]
                seekposition += 3*nbytes_Float64
            end
        end
        seek(fp,seekposition)
        index = read(fp,Int64)
        numatoms = read(fp,Int64)
        @assert numatoms == trjfile.numatoms_list[fileindex] "data structure might be wrong"
        @assert index == itrj
        #println("\t $index $numatoms")

        structure = Array{Atom,1}(undef,numatoms)

        #structure = zeros(Float64,3,numatoms)
        for iatom=1:numatoms
            x = read(fp,Float64)
            y = read(fp,Float64)
            z = read(fp,Float64)
            structure[iatom] = Atom(Float64[x,y,z])
        end
        seek(fp,0)
        close(fp)
        #println(fileindex)

        positions = get_relative_positions(Snapshot(structure))

        for iatom=1:numatoms
            structure[iatom] = Atom(positions[:,iatom])
        end

        return Snapshot(structure),fileindex

    end


    function calc_MSD_fromfile(inputdirname,dirname,itrj_start,
                        period,start_itrj,maxsteps)

        trjfile = Trjfile(inputdirname)
        println(trjfile.numdata_list)
        println(trjfile.index_list)
        println(trjfile.numatoms_list)



        snapshot,fileindex = get_structure(trjfile,itrj_start)
        snapshot_old = deepcopy(snapshot)

        numatoms = trjfile.numatoms_list[fileindex]

        count = 0

        positions = zeros(Float64,3,numatoms,maxsteps)
        positions_old = zeros(Float64,3,numatoms)

        lastindex = get_lastindex(trjfile)

        msd = zeros(maxsteps)
        phi = zeros(numatoms,maxsteps)
        ri = zeros(3)
        drj = zeros(3)
        dphi = 0


        count = 0
        for itrj=itrj_start:period:lastindex

            if itrj+maxsteps-1 > lastindex
                break
            end
            count += 1
            snapshot,fileindex = get_structure(trjfile,itrj)

            for it=2:maxsteps
                jtrj = itrj+it-1
                snapshot_j,fileindex_j = get_structure(trjfile,jtrj)

                for iatom=1:numatoms
                    dphi_x = (snapshot_j[iatom].x*a0-snapshot[iatom].x*a0)^2
                    dphi_y = (snapshot_j[iatom].y*a0-snapshot[iatom].y*a0)^2
                    dphi_z = (snapshot_j[iatom].z*a0-snapshot[iatom].z*a0)^2
                    dphi = dphi_x+dphi_y+dphi_z
                    msd[it] += dphi
                    phi[iatom,it] += dphi
                end
            end

            #println(msd[2],"\t",count)

            fp = open(dirname*"/atomic_MSD_from$(itrj_start)-th_file.txt","w")
            fp2 = open(dirname*"/MSD_from$(itrj_start)-th_file.txt","w")


            for it = 1:maxsteps
                println(fp2,it,"\t",msd[it]/(count*numatoms))
                print(fp,it,"\t")
                for iatom=1:numatoms
                    print(fp,phi[iatom,it]/count,"\t")
                end
                println(fp,"\t")
            end
            close(fp)
    
            close(fp2)
            println("itrj = ",itrj)


        end



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
            if itrj+maxsteps-1 > length(trj)
                break
            end
            count += 1

            for iatom=1:trj.numatoms
                for i=1:3
                    ri[i] = positions[i,iatom,itrj]*a0
                end
                for it = 1:maxsteps
                    jtrj = itrj+it-1
                    
                    dphi = 0.0
                    for i=1:3
                        dphi += (positions[i,iatom,jtrj]*a0 - ri[i])^2
                        #drj[i] = positions[i,iatom,jtrj]*a0 - ri[i]
                    end
                    msd[it] += dphi
                    phi[iatom,it] += dphi
                end                

            end
            #println(msd[2],"\t",count)
            fp = open(dirname*"/atomic_MSD_from$(start_itrj)-th.txt","w")
            fp2 = open(dirname*"/MSD_from$(start_itrj)-th.txt","w")


            for it = 1:maxsteps

                #println(positions[:,:,istart+1])
                #return

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
