using Revise
using ObsIO
using BenchmarkTools

ipath = "/home/work/b-physics/data/J307/J307r000.in"
dpath = "/home/work/b-physics/data/J307/b1/J307r000_cnfg1.mesons.dat"
corr =  ObsIO.read_input_file(ipath,ObsIO.G5,ObsIO.G5)
cdata = ObsIO.read_mesons(dpath,ObsIO.G5,ObsIO.G5)
corr = ObsIO.corr_obs.(cdata,corr,L=1,real=true)
write_corr("test.bdio",corr[1],override=true)
C = ObsIO._read_corr("test.bdio")

#Periodic boundary conditions
ipath = "/home/work/b-physics/data/B450/B450r000.in"
corr = @benchmark ObsIO.read_input_file(ipath,ObsIO.G5,ObsIO.G5)


0==0
