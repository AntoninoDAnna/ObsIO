using Revise
using ObsIO

path = joinpath(homedir(),"ift_desk/codes/input_parser/H101DeltaT.in")

Pr,corr = ObsIO.read_input_file(path)
