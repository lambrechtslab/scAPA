#Create a softlink for scAPA program.
run(`ln -s $(@__DIR__)/bin/scAPA $(ENV["HOME"]*"/")`)
