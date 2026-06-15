################################################################################
## EXTRACT PATHWAY ANNOTATIONS FROM YEAST-GEM MODELS
################################################################################

using SBML, DataFrames, CSV

# Load S. cerevisiae model (reference).
# Download Saccharomyces_cerevisiae.xml from the Yeast-Species-GEMs dataset
# (Shen et al. 2018, https://github.com/SysBioChalmers/Yeast-Species-GEMs)
# and place it at data/Saccharomyces_cerevisiae.xml, or set the path below.
REPO_ROOT  = joinpath(@__DIR__, "..", "..")
model_path = get(ENV, "YEAST_GEM_XML", joinpath(REPO_ROOT, "data", "Saccharomyces_cerevisiae.xml"))

println("Loading model: ", model_path)
model = readSBML(model_path)

# Extract metabolites with their IDs
metabolites = DataFrame(
    metabolite_id = String[],
    name = String[]
)

for (id, met) in model.species
    push!(metabolites, (id, met.name))
end

println("\nExtracted ", nrow(metabolites), " metabolites")

# Extract reactions with subsystem annotations
reactions = DataFrame(
    reaction_id = String[],
    name = String[],
    subsystem = String[]
)

for (id, rxn) in model.reactions
    subsystem = haskey(rxn.notes, "SUBSYSTEM") ? rxn.notes["SUBSYSTEM"] : "Unknown"
    push!(reactions, (id, rxn.name, subsystem))
end

println("Extracted ", nrow(reactions), " reactions")

# Count metabolites per subsystem (via reactions they participate in)
metabolite_subsystems = DataFrame(
    metabolite_id = String[],
    subsystem = String[],
    n_reactions = Int[]
)

for (rxn_id, rxn) in model.reactions
    subsystem = haskey(rxn.notes, "SUBSYSTEM") ? rxn.notes["SUBSYSTEM"] : "Unknown"

    # Get metabolites from this reaction (reactants + products)
    met_ids = Set{String}()
    for (met_id, stoich) in rxn.reactants
        push!(met_ids, met_id)
    end
    for (met_id, stoich) in rxn.products
        push!(met_ids, met_id)
    end

    # Record metabolite-subsystem associations
    for met_id in met_ids
        push!(metabolite_subsystems, (met_id, subsystem, 1))
    end
end

# Aggregate: count reactions per metabolite-subsystem pair
met_sub_summary = combine(groupby(metabolite_subsystems, [:metabolite_id, :subsystem]),
                          :n_reactions => sum => :n_reactions)

println("\nMetabolite-subsystem associations: ", nrow(met_sub_summary))

# Save outputs
output_dir = joinpath(REPO_ROOT, "results", "acr")
CSV.write(joinpath(output_dir, "metabolite_pathway_annotations.csv"), met_sub_summary)
CSV.write(joinpath(output_dir, "reaction_subsystems.csv"), reactions)

println("\nSaved pathway annotations to:")
println("  ", joinpath(output_dir, "metabolite_pathway_annotations.csv"))
println("  ", joinpath(output_dir, "reaction_subsystems.csv"))

# Preview subsystems
println("\n\nUnique subsystems:")
unique_subs = sort(unique(met_sub_summary.subsystem))
for (i, sub) in enumerate(unique_subs)
    println("  ", i, ". ", sub)
end
