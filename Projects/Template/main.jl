push!(LOAD_PATH, pwd()*"/Modules/")
#push!(LOAD_PATH, pwd()*"/Projects/ProjectX/Local Modules/") # Uncomment this line if you want to use the local modules for project X

using Revise, Plots, PlotlyJS # Always leave this on!
Plots.plotlyjs() # Set PlotlyJS as the backend for Plots.jl

#### Other modules that you might need go here ####

###################################################

#### Your code goes here ####