using Revise
using ObsIO

ipath = joinpath(homedir(),"ift_desk/data/J307/J307r000.in")
dpath = joinpath(homedir(),"ift_desk/data/J307/b1/J307r000_cnfg1.mesons.dat")
corr = ObsIO.read_input_file(ipath,ObsIO.G5,ObsIO.G5)
cdata = ObsIO.read_mesons(dpath,ObsIO.G5,ObsIO.G5)
corr = ObsIO.corr_obs.(cdata,corr,L=1,real=true)
