"""
Giant Kinetic Module Reaction Presence Analysis - OPTIMIZED VERSION

Key optimizations:
1. Faster XML parsing for pathways (streaming regex, no full model load)
2. Better progress reporting with timing information
3. Optimized memory usage

Note: Incidence matrices cannot be cached since random_0 models change between runs.

Usage:
    julia --project=@analysis --threads=auto giant_module_reaction_analysis_optimized.jl
"""

using JLD2, CSV, DataFrames, Statistics
using SBMLFBCModels, AbstractFBCModels, COBREXA
using COCOA
using SparseArrays
using ProgressMeter

# =============================================================================
# Configuration
# =============================================================================

# Paths relative to repo root.
# JLD2 kinetic results: extract random_0.tar.gz from Zenodo into data/jld2/random_0/
# PRPD models: extract prpd_models/random_0.tar.gz from Zenodo into data/prpd_models/random_0/
const REPO_ROOT       = joinpath(@__DIR__, "..", "..")
const KINETIC_DIR     = get(ENV, "JLD2_RESULTS_DIR", joinpath(REPO_ROOT, "data", "jld2", "random_0"))
const MODELS_DIR      = get(ENV, "PRPD_MODELS_DIR", joinpath(REPO_ROOT, "data", "prpd_models", "random_0"))
const OUTPUT_DIR      = joinpath(REPO_ROOT, "results", "kinetic")

mkpath(OUTPUT_DIR)

# =============================================================================
# Utility Functions
# =============================================================================

"""
    extract_original_reaction(rxn_id::Symbol) -> String

Extract the original reaction ID from an elementary step reaction ID.
"""
function extract_original_reaction(rxn_id::Symbol)
    id_str = String(rxn_id)
    m = match(r"^(.+)_E\d+_(SB\d+|CAT|PR\d+)$", id_str)
    if m !== nothing
        return m.captures[1]
    else
        return id_str
    end
end

"""
    extract_step_type(rxn_id::Symbol) -> String

Extract the elementary step type (SB, CAT, PR, or "original" for non-decomposed).
"""
function extract_step_type(rxn_id::Symbol)
    id_str = String(rxn_id)
    m = match(r"_E\d+_(SB|CAT|PR)\d*$", id_str)
    if m !== nothing
        return m.captures[1]
    else
        return "original"
    end
end

"""
    extract_pathways_from_sbml_fast(sbml_path::String) -> Dict{String, Vector{String}}

Fast pathway extraction using streaming regex (no full XML parsing).
"""
function extract_pathways_from_sbml_fast(sbml_path::String)
    reaction_to_pathways = Dict{String,Vector{String}}()

    # Read file in chunks to handle large files
    content = read(sbml_path, String)

    # Find groups section (much faster than full XML parsing)
    groups_match = match(r"<groups:listOfGroups>(.*?)</groups:listOfGroups>"s, content)
    if groups_match === nothing
        return reaction_to_pathways
    end

    groups_section = groups_match.captures[1]

    # Extract each group
    group_pattern = r"<groups:group[^>]+groups:name=\"([^\"]+)\"[^>]*>(.*?)</groups:group>"s

    for group_match in eachmatch(group_pattern, groups_section)
        pathway_name = group_match.captures[1]
        group_content = group_match.captures[2]

        # Find all member references
        member_pattern = r"groups:idRef=\"([^\"]+)\""

        for member_match in eachmatch(member_pattern, group_content)
            reaction_id = member_match.captures[1]

            if !haskey(reaction_to_pathways, reaction_id)
                reaction_to_pathways[reaction_id] = String[]
            end
            push!(reaction_to_pathways[reaction_id], pathway_name)
        end
    end

    return reaction_to_pathways
end

# =============================================================================
# Core Processing
# =============================================================================

"""
    process_species(species_id::String) -> NamedTuple

Process a single species using cached incidence matrices.
"""
function process_species(species_id::String)
    # 1. Load kinetic module results
    kinetic_path = joinpath(KINETIC_DIR, "kinetic_modules_$(species_id).jld2")
    data = JLD2.load(kinetic_path)

    kinetic_modules = get(data, "kinetic_modules", Vector{Set{Symbol}}())
    stats = get(data, "statistics", Dict())
    n_total_complexes = get(stats, "n_total_complexes", 0)

    if isempty(kinetic_modules)
        return (
            species_id=species_id,
            giant_reactions=Set{Symbol}(),
            reaction_to_original=Dict{Symbol,String}(),
            giant_module_size=0,
            n_total_complexes=n_total_complexes
        )
    end

    # 2. Find giant module
    module_sizes = length.(kinetic_modules)
    giant_idx = argmax(module_sizes)
    giant_module = kinetic_modules[giant_idx]
    giant_module_size = module_sizes[giant_idx]

    # 3. Load model and compute incidence matrix
    model_path = joinpath(MODELS_DIR, "$(species_id).xml")
    model = COBREXA.load_model(model_path)
    A, complex_ids, reaction_ids = COCOA.incidence(model; return_ids=true)

    # 4. Build complex index lookup
    complex_to_idx = Dict(id => i for (i, id) in enumerate(complex_ids))

    # 5. Find substrate reactions
    giant_complex_indices = Int[]
    n_missing = 0

    for cplx_id in giant_module
        idx = get(complex_to_idx, cplx_id, 0)
        if idx == 0
            n_missing += 1
        else
            push!(giant_complex_indices, idx)
        end
    end

    if n_missing > 0
        @warn "$(species_id): $(n_missing)/$(giant_module_size) giant module complexes not found"
    end

    # 6. Extract substrate reactions
    giant_reactions = Set{Symbol}()
    for idx in giant_complex_indices
        row = A[idx, :]
        for (j, v) in zip(findnz(row)...)
            if v == -1
                push!(giant_reactions, reaction_ids[j])
            end
        end
    end

    # 7. Build reaction-to-original mapping
    reaction_to_original = Dict{Symbol,String}(
        rxn_id => extract_original_reaction(rxn_id) for rxn_id in giant_reactions
    )

    return (
        species_id=species_id,
        giant_reactions=giant_reactions,
        reaction_to_original=reaction_to_original,
        giant_module_size=giant_module_size,
        n_total_complexes=n_total_complexes
    )
end

# =============================================================================
# Matrix Construction
# =============================================================================

"""
    build_reaction_presence_matrix(species_results, species_order) -> DataFrame
"""
function build_reaction_presence_matrix(
    species_results::Dict{String,Set{Symbol}},
    species_order::Vector{String}
)
    all_reactions = Set{Symbol}()
    for rxns in values(species_results)
        union!(all_reactions, rxns)
    end
    all_reactions_sorted = sort(collect(all_reactions))
    reaction_to_row = Dict(rxn => i for (i, rxn) in enumerate(all_reactions_sorted))

    println("  Found $(length(all_reactions_sorted)) unique reactions across $(length(species_order)) species")

    n_reactions = length(all_reactions_sorted)
    n_species = length(species_order)
    matrix = zeros(Int, n_reactions, n_species)

    for (j, species) in enumerate(species_order)
        if haskey(species_results, species)
            for rxn in species_results[species]
                i = reaction_to_row[rxn]
                matrix[i, j] = 1
            end
        end
    end

    df = DataFrame(matrix, Symbol.(species_order))
    insertcols!(df, 1, :reaction_id => String.(all_reactions_sorted))

    return df
end

# =============================================================================
# Main
# =============================================================================

function main()
    println("="^70)
    println("Giant Kinetic Module Reaction Presence Analysis - OPTIMIZED")
    println("="^70)

    # Discover species
    jld2_files = filter(f -> endswith(f, ".jld2"), readdir(KINETIC_DIR))
    species_list = [replace(f, "kinetic_modules_" => "", ".jld2" => "") for f in jld2_files]

    species_list = filter(species_list) do sp
        isfile(joinpath(MODELS_DIR, "$(sp).xml"))
    end
    sort!(species_list)

    println("\nFound $(length(species_list)) species")
    println("Using $(Threads.nthreads()) threads")

    # Process all species
    println("\n1. Processing species (loading kinetic results + incidence matrices)...")
    results = Vector{Union{Nothing,NamedTuple}}(undef, length(species_list))

    p = Progress(length(species_list); desc="Processing: ", dt=0.5)
    Threads.@threads for i in eachindex(species_list)
        results[i] = try
            r = process_species(species_list[i])
            next!(p)
            r
        catch e
            @warn "Failed for $(species_list[i]): $e"
            next!(p)
            nothing
        end
    end
    finish!(p)

    valid_results = filter(!isnothing, results)
    println("  Successfully processed $(length(valid_results))/$(length(species_list)) species")

    # Build data structures
    species_reactions = Dict{String,Set{Symbol}}()
    all_reaction_to_original = Dict{Symbol,String}()
    summary_rows = NamedTuple[]

    for r in valid_results
        species_reactions[r.species_id] = r.giant_reactions
        merge!(all_reaction_to_original, r.reaction_to_original)

        n_original = length(unique(values(r.reaction_to_original)))
        push!(summary_rows, (
            species_id=r.species_id,
            n_giant_reactions=length(r.giant_reactions),
            n_original_reactions=n_original,
            giant_module_size=r.giant_module_size,
            n_total_complexes=r.n_total_complexes,
            giant_proportion=r.n_total_complexes > 0 ? r.giant_module_size / r.n_total_complexes : 0.0
        ))
    end

    species_order = sort(collect(keys(species_reactions)))

    # 2. Build presence matrix
    println("\n2. Building reaction presence/absence matrix...")
    presence_matrix = build_reaction_presence_matrix(species_reactions, species_order)

    # 3. Extract pathway information
    println("\n3. Extracting pathway information from original models...")
    species_pathways = Dict{String,Dict{String,Vector{String}}}()

    pathway_results = Vector{Union{Nothing,Tuple{String,Dict{String,Vector{String}}}}}(undef, length(species_order))

    p2 = Progress(length(species_order); desc="Loading pathways: ", dt=0.5)
    Threads.@threads for i in eachindex(species_order)
        species = species_order[i]
        original_model_path = joinpath(ORIGINAL_GEMS_DIR, "$(species).xml")

        pathway_results[i] = if isfile(original_model_path)
            try
                pathways = extract_pathways_from_sbml_fast(original_model_path)
                next!(p2)
                (species, pathways)
            catch e
                @warn "Failed to extract pathways for $(species): $e"
                next!(p2)
                (species, Dict{String,Vector{String}}())
            end
        else
            @warn "Original model not found for $(species)"
            next!(p2)
            (species, Dict{String,Vector{String}}())
        end
    end
    finish!(p2)

    for result in pathway_results
        if result !== nothing
            species, pathways = result
            species_pathways[species] = pathways
        end
    end

    println("  Extracted pathways for $(count(sp -> !isempty(species_pathways[sp]), species_order)) species")

    # 4. Build reaction metadata
    println("\n4. Building reaction metadata...")
    n_species = length(species_order)
    metadata_rows = Vector{NamedTuple}(undef, nrow(presence_matrix))

    species_cols = [Symbol(sp) for sp in species_order]

    Threads.@threads for idx in 1:nrow(presence_matrix)
        row = presence_matrix[idx, :]
        rxn_id = row[:reaction_id]
        rxn_sym = Symbol(rxn_id)

        n_present = 0
        for col in species_cols
            n_present += row[col]
        end

        original_rxn_id = get(all_reaction_to_original, rxn_sym, rxn_id)

        pathways_set = Set{String}()
        for (i, species) in enumerate(species_order)
            if row[species_cols[i]] == 1
                if haskey(species_pathways, species)
                    rxn_pathways = get(species_pathways[species], original_rxn_id, String[])
                    union!(pathways_set, rxn_pathways)
                end
            end
        end
        pathways_str = isempty(pathways_set) ? "" : join(sort(collect(pathways_set)), "; ")

        metadata_rows[idx] = (
            reaction_id=rxn_id,
            original_reaction_id=original_rxn_id,
            step_type=extract_step_type(rxn_sym),
            pathways=pathways_str,
            n_species_present=n_present,
            proportion_present=n_present / n_species
        )
    end

    metadata_df = DataFrame(metadata_rows)
    sort!(metadata_df, :n_species_present, rev=true)

    # 5. Build species summary
    println("\n5. Building species summary...")
    summary_df = DataFrame(summary_rows)
    sort!(summary_df, :species_id)

    # 6. Save outputs
    println("\n6. Saving outputs...")

    matrix_path = joinpath(OUTPUT_DIR, "giant_module_reaction_presence_matrix.csv")
    CSV.write(matrix_path, presence_matrix)
    println("  Saved: $matrix_path ($(nrow(presence_matrix)) reactions × $(ncol(presence_matrix)-1) species)")

    metadata_path = joinpath(OUTPUT_DIR, "giant_module_reaction_metadata.csv")
    CSV.write(metadata_path, metadata_df)
    println("  Saved: $metadata_path")

    summary_path = joinpath(OUTPUT_DIR, "giant_module_reaction_summary.csv")
    CSV.write(summary_path, summary_df)
    println("  Saved: $summary_path")

    # 7. Print summary statistics
    println("\n" * "="^70)
    println("Summary Statistics")
    println("="^70)
    println("  Species processed: $(nrow(summary_df))")
    println("  Unique reactions in any giant module: $(nrow(presence_matrix))")
    println("  Unique original reactions: $(length(unique(metadata_df.original_reaction_id)))")
    println("\n  Giant module reactions per species:")
    println("    Mean:   $(round(mean(summary_df.n_giant_reactions), digits=1))")
    println("    Median: $(round(median(summary_df.n_giant_reactions), digits=1))")
    println("    Min:    $(minimum(summary_df.n_giant_reactions))")
    println("    Max:    $(maximum(summary_df.n_giant_reactions))")
    println("\n  Reaction conservation:")
    println("    Present in >90% species: $(sum(metadata_df.proportion_present .> 0.9))")
    println("    Present in >50% species: $(sum(metadata_df.proportion_present .> 0.5))")
    println("    Present in >15% species: $(sum(metadata_df.proportion_present .> 0.15))")
    println("    Present in <5% species:  $(sum(metadata_df.proportion_present .< 0.05))")

    println("\n  Step type breakdown:")
    for g in groupby(metadata_df, :step_type)
        println("    $(g.step_type[1]): $(nrow(g)) reactions")
    end

    println("\n  Pathway annotation:")
    n_with_pathways = sum(!isempty(p) for p in metadata_df.pathways)
    println("    Reactions with pathway annotation: $(n_with_pathways)/$(nrow(metadata_df)) ($(round(100*n_with_pathways/nrow(metadata_df), digits=1))%)")

    all_pathways = Set{String}()
    for pathways_str in metadata_df.pathways
        if !isempty(pathways_str)
            for pathway in split(pathways_str, "; ")
                push!(all_pathways, pathway)
            end
        end
    end
    println("    Unique pathways represented: $(length(all_pathways))")

    println("\nDone!")
end

main()
