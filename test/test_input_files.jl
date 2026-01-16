using Revise
using ObsIO

path = "/home/antonino/ift_desk/codes/input_parser/H101DeltaT.in"

props,points = ObsIO.read_input_file(path)
