using PIMDtools
using Test

function test1()
    itrj_start = 4
    itrj_end = 100
    atominfo = [("Pd",20),("Ru",16),("Al",92)]
    a = 12.6238002777
    unitcell = Array{Array{Float64,1},1}(undef,3)
    for i=1:3
        unitcell[i] = zeros(Float64,3)
        unitcell[i][i] = a/a0
    end
    trj = Trajectory("trj.out",itrj_start,itrj_end,atominfo,unitcell=unitcell)
    positions = get_relative_positions(trj)

    dirname = "./"
    period = 10
    start_itrj = itrj_start 
    maxsteps = 10
    calc_MSD(trj,dirname,positions,period,start_itrj,maxsteps)
end

function testbinary()
    itrj_start = 4
    itrj_end = 100
    atominfo = [("Pd",20),("Ru",16),("Al",92)]
    a = 12.6238002777
    mkpath("testdir")
    filename = "trj.out"
    dirname = "./testdir"
    num_of_structures = 5
    

    inputdirname = "./testdir"
    period = 10
    start_itrj = itrj_start 
    maxsteps = 10

    unitcell = Array{Array{Float64,1},1}(undef,3)
    for i=1:3
        unitcell[i] = zeros(Float64,3)
        unitcell[i][i] = a/a0
    end

    toBinary_r(filename,dirname,itrj_start,itrj_end,num_of_structures,atominfo,unitcell)


    calc_MSD_fromfile(inputdirname,"./",itrj_start,
        period,start_itrj,maxsteps)
end



@testset "PIMDtools.jl" begin
    
    test1()
    testbinary()
    # Write your tests here.
end


