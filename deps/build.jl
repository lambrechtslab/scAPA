#Create a softlink for scAPA program.
fn=joinpath(ENV["HOME"], "scAPA")
if isfile(fn)
    run(`unlink $fn`)
end
run(`ln -s $(@__DIR__)/../bin/scAPA $fn`)
