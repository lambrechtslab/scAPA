    parseUserInputs()
    
    ## Extract junctions
    fetch_junction_reads()
    union_junctions()
    add_unique_junction_id()

    ## Construct ASEs
    ase_info_file=joinpath(OUTDIR, "valid_splicingEvents_allAsTyp.jld2")