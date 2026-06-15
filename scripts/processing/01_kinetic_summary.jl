"""
Kinetic Module & ACR Phylogenetic Analysis

Replicates key analyses from Langary et al. 2025 and extends with phylogenetic analysis.
- Giant module size comparison across yeast species
- ACR/ACRR counts per species
- Phylogenetic signal in ACR metabolite sets
- Ancestral state reconstruction for common ACR metabolites

Usage:
    julia --threads=auto kinetic_phylogenetic_analysis.jl
"""

using JLD2, CSV, DataFrames, Statistics
using CairoMakie
using ProgressMeter
using Distances

# =============================================================================
# Configuration
# =============================================================================

# Paths relative to repo root. JLD2 kinetic results must be extracted from
# the Zenodo archive (random_0.tar.gz) into a local directory first.
const REPO_ROOT     = joinpath(@__DIR__, "..", "..")
const RESULTS_DIR   = get(ENV, "JLD2_RESULTS_DIR", joinpath(REPO_ROOT, "data", "jld2", "random_0"))
const METADATA_FILE = joinpath(REPO_ROOT, "data", "annotations", "cocoa_results_merged.csv")
const TREE_FILE     = joinpath(REPO_ROOT, "data", "phylogeny", "332_2408OGs_timetree_mcmctree.nwk")
const OUTPUT_DIR    = joinpath(REPO_ROOT, "results", "kinetic")

# =============================================================================
# Metabolite ID Utilities
# =============================================================================

"""
    decode_id(encoded_id::String) -> String

Decode special characters in metabolite IDs.
Common encodings:
- __91__ = [
- __93__ = ]
- __45__ = -
- __43__ = +
"""
function decode_id(encoded_id::String)
    decoded = string(encoded_id)
    decoded = replace(decoded, "__91__" => "[")
    decoded = replace(decoded, "__93__" => "]")
    decoded = replace(decoded, "__45__" => "-")
    decoded = replace(decoded, "__43__" => "+")
    return decoded
end

"""
    is_real_metabolite(met_id) -> Bool

Check if a metabolite ID represents a real metabolite (not an enzyme complex
or intermediate from ordered/random binding expansion).

Filters out:
- CPLX_* : enzyme-substrate complexes
- E\\d+   : enzyme identifiers (e.g., E354, E479)
- *+*    : bound states (contain + in ID)
"""
function is_real_metabolite(met_id)
    id_str = string(met_id)

    # Filter out enzyme complexes (CPLX_)
    if occursin("CPLX", id_str)
        return false
    end

    # Filter out enzyme identifiers (E followed by digits)
    if occursin(r"_E\d+", id_str) || startswith(id_str, "E") && occursin(r"^E\d+", id_str)
        return false
    end

    # Filter out bound states (contain + which indicates enzyme-metabolite binding)
    if occursin("+", id_str) || occursin("__43__", id_str)
        return false
    end

    return true
end

"""
    filter_real_metabolites(metabolites::Vector) -> Vector

Filter a vector of metabolite IDs to keep only real metabolites.
"""
function filter_real_metabolites(metabolites::Vector)
    return filter(is_real_metabolite, metabolites)
end

# =============================================================================
# Data Loading Functions
# =============================================================================

"""
    load_species_metadata(csv_path) -> DataFrame

Load species metadata with clade assignments.
"""
function load_species_metadata(csv_path::String)
    df = CSV.read(csv_path, DataFrame)
    # Key columns: old_species_id (col 49), Major clade (col 51)
    select!(df,
        "old_species_id" => :species_id,
        "Species name" => :species_name,
        "Major clade" => :clade,
        "n_complexes" => :n_complexes_metadata
    )
    return df
end

"""
    load_kinetic_result(filepath) -> NamedTuple

Load a single kinetic results JLD2 file and extract key metrics.
"""
function load_kinetic_result(filepath::String)
    data = JLD2.load(filepath)

    # Extract species name from filename
    basename_str = basename(filepath)
    species_id = replace(basename_str, "kinetic_modules_" => "", ".jld2" => "")

    # Get statistics dict
    stats = get(data, "statistics", Dict())

    # Get ACR metabolites and ACRR pairs (raw, including enzyme complexes)
    acr_mets_raw = get(data, "acr_metabolites", Symbol[])
    acrr_pairs_raw = get(data, "acrr_pairs", Tuple{Symbol,Symbol}[])
    kinetic_modules = get(data, "kinetic_modules", Vector{Set{Symbol}}())

    # Filter to keep only real metabolites (exclude CPLX, E###, bound states)
    acr_mets = filter_real_metabolites(acr_mets_raw)
    acrr_pairs = filter(p -> is_real_metabolite(p[1]) && is_real_metabolite(p[2]), acrr_pairs_raw)

    # Calculate metrics (using filtered counts)
    n_acr = length(acr_mets)
    n_acrr = length(acrr_pairs)
    n_acr_raw = length(acr_mets_raw)  # Keep raw count for reference
    largest_module_size = get(stats, "largest_module_size", 0)
    n_total_complexes = get(stats, "n_total_complexes", 0)
    n_kinetic_modules = get(stats, "n_kinetic_modules", length(kinetic_modules))

    # Proportion of elementary reaction steps with substrate complex in the giant module
    # (= n_giant_reactions / n_total_split_reactions, computed by joining with
    # giant_module_reaction_summary.csv.gz and split_model_sizes.csv.gz).
    # The value stored here is complex-based (giant_module_size / n_total_complexes)
    # and is overwritten by fix_giant_proportion.py before committing to the repo.
    proportion = n_total_complexes > 0 ? largest_module_size / n_total_complexes : 0.0

    return (
        species_id = species_id,
        n_acr = n_acr,
        n_acrr = n_acrr,
        n_acr_raw = n_acr_raw,
        largest_module_size = largest_module_size,
        n_total_complexes = n_total_complexes,
        n_kinetic_modules = n_kinetic_modules,
        giant_proportion = proportion,
        acr_metabolites = acr_mets,
        acrr_pairs = acrr_pairs
    )
end

"""
    load_all_kinetic_results(results_dir) -> DataFrame, Dict

Load all kinetic results and return summary DataFrame plus ACR metabolites dict.
"""
function load_all_kinetic_results(results_dir::String)
    jld2_files = filter(f -> endswith(f, ".jld2"), readdir(results_dir, join=true))

    println("Loading $(length(jld2_files)) kinetic result files...")

    results = []
    acr_dict = Dict{String, Vector{Symbol}}()  # species_id => ACR metabolites

    @showprogress for filepath in jld2_files
        try
            result = load_kinetic_result(filepath)
            push!(results, (
                species_id = result.species_id,
                n_acr = result.n_acr,
                n_acrr = result.n_acrr,
                n_acr_raw = result.n_acr_raw,  # Include raw count (with CPLX/E###)
                largest_module_size = result.largest_module_size,
                n_total_complexes = result.n_total_complexes,
                n_kinetic_modules = result.n_kinetic_modules,
                giant_proportion = result.giant_proportion
            ))
            acr_dict[result.species_id] = result.acr_metabolites
        catch e
            @warn "Failed to load $filepath: $e"
        end
    end

    df = DataFrame(results)
    println("Loaded $(nrow(df)) species successfully")
    println("   Note: ACR counts filtered to exclude CPLX/enzyme complexes")
    println("   Raw ACR range: $(minimum(df.n_acr_raw))-$(maximum(df.n_acr_raw)), Filtered: $(minimum(df.n_acr))-$(maximum(df.n_acr))")

    return df, acr_dict
end

# =============================================================================
# ACR Matrix Construction
# =============================================================================

"""
    build_acr_presence_matrix(acr_dict, species_order) -> DataFrame

Build presence/absence matrix of ACR metabolites across species.
Rows = metabolites, Columns = species
"""
function build_acr_presence_matrix(acr_dict::Dict{String, Vector{Symbol}}, species_order::Vector{String})
    # Get all unique ACR metabolites across all species
    all_mets = Set{Symbol}()
    for mets in values(acr_dict)
        union!(all_mets, mets)
    end
    all_mets = sort(collect(all_mets))

    println("Found $(length(all_mets)) unique ACR metabolites across $(length(species_order)) species")

    # Build matrix
    matrix = zeros(Int, length(all_mets), length(species_order))

    for (j, species) in enumerate(species_order)
        if haskey(acr_dict, species)
            for met in acr_dict[species]
                i = findfirst(==(met), all_mets)
                if i !== nothing
                    matrix[i, j] = 1
                end
            end
        end
    end

    # Convert to DataFrame
    df = DataFrame(matrix, Symbol.(species_order))
    insertcols!(df, 1, :metabolite => String.(all_mets))

    return df
end

"""
    get_metabolite_conservation(acr_matrix) -> DataFrame

Calculate conservation statistics for each ACR metabolite.
"""
function get_metabolite_conservation(acr_matrix::DataFrame)
    n_species = ncol(acr_matrix) - 1  # Exclude metabolite column

    conservation = DataFrame(
        metabolite = acr_matrix.metabolite,
        n_species_acr = [sum(row[2:end]) for row in eachrow(acr_matrix)],
        proportion = [sum(row[2:end]) / n_species for row in eachrow(acr_matrix)]
    )

    sort!(conservation, :n_species_acr, rev=true)
    return conservation
end

# =============================================================================
# Phylogenetic Analysis
# =============================================================================

"""
    calculate_jaccard_similarity(set1, set2) -> Float64

Calculate Jaccard similarity between two sets.
"""
function jaccard_similarity(set1::Set, set2::Set)
    if isempty(set1) && isempty(set2)
        return 1.0  # Both empty = identical
    end
    intersection = length(intersect(set1, set2))
    union_size = length(union(set1, set2))
    return intersection / union_size
end

"""
    calculate_all_jaccard_similarities(acr_dict, species_order) -> Matrix

Calculate pairwise Jaccard similarities between all species ACR sets.
"""
function calculate_all_jaccard_similarities(acr_dict::Dict{String, Vector{Symbol}}, species_order::Vector{String})
    n = length(species_order)
    sim_matrix = zeros(n, n)

    for i in 1:n
        set_i = Set(get(acr_dict, species_order[i], Symbol[]))
        for j in i:n
            set_j = Set(get(acr_dict, species_order[j], Symbol[]))
            sim = jaccard_similarity(set_i, set_j)
            sim_matrix[i, j] = sim
            sim_matrix[j, i] = sim
        end
    end

    return sim_matrix
end

# =============================================================================
# Simple Newick Parser (for basic tree operations)
# =============================================================================

"""
Simple Newick tree parser - extracts tip names and basic structure.
For full phylogenetic operations, use Phylo.jl
"""
function parse_newick_tips(newick_str::String)
    # Extract tip names (labels before : or , or ) )
    tips = String[]

    # Remove comments in brackets [&...]
    clean = replace(newick_str, r"\[&[^\]]*\]" => "")

    # Find all labels
    for m in eachmatch(r"([A-Za-z0-9_]+):", clean)
        push!(tips, m.captures[1])
    end

    return tips
end

# =============================================================================
# Plotting Functions (CairoMakie)
# =============================================================================

"""
    plot_giant_module_by_clade(df) -> Figure

Box/violin plot of giant module proportion by clade.
"""
function plot_giant_module_by_clade(df::DataFrame; output_path::String="")
    # Get unique clades and sort by median proportion
    clade_medians = combine(groupby(df, :clade), :giant_proportion => median => :median_prop)
    sort!(clade_medians, :median_prop, rev=true)
    clade_order = clade_medians.clade

    # Create figure
    fig = Figure(size=(1200, 800))
    ax = Axis(fig[1, 1],
        xlabel = "Giant module proportion (fraction of complexes)",
        ylabel = "Clade",
        title = "Giant Kinetic Module Size by Clade",
        yticks = (1:length(clade_order), clade_order)
    )

    # Plot boxplots
    for (i, clade) in enumerate(clade_order)
        clade_data = filter(row -> row.clade == clade, df).giant_proportion
        if !isempty(clade_data)
            boxplot!(ax, fill(i, length(clade_data)), clade_data,
                orientation=:horizontal, color=Makie.wong_colors()[mod1(i, 7)])
        end
    end

    if !isempty(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

"""
    plot_acr_counts_by_species(df; sort_by_clade=true) -> Figure

Horizontal bar plot of ACR metabolite counts per species.
"""
function plot_acr_counts_by_species(df::DataFrame; output_path::String="", max_species::Int=50)
    # Sort by clade then by n_acr
    df_sorted = sort(df, [:clade, :n_acr], rev=[false, true])

    # Limit to top species if too many
    if nrow(df_sorted) > max_species
        # Take top species from each clade
        df_sorted = combine(groupby(df_sorted, :clade)) do sdf
            first(sort(sdf, :n_acr, rev=true), min(5, nrow(sdf)))
        end
    end

    # Create color mapping for clades
    clades = unique(df_sorted.clade)
    clade_colors = Dict(c => Makie.wong_colors()[mod1(i, 7)] for (i, c) in enumerate(clades))

    fig = Figure(size=(1000, max(600, nrow(df_sorted) * 15)))
    ax = Axis(fig[1, 1],
        xlabel = "Number of ACR metabolites",
        ylabel = "Species",
        title = "Absolute Concentration Robustness by Species",
        yticks = (1:nrow(df_sorted), df_sorted.species_id),
        yticklabelsize = 8
    )

    colors = [clade_colors[c] for c in df_sorted.clade]
    barplot!(ax, 1:nrow(df_sorted), df_sorted.n_acr,
        direction=:x, color=colors)

    if !isempty(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

"""
    plot_acrr_counts_by_species(df) -> Figure

Horizontal bar plot of ACRR pair counts per species.
"""
function plot_acrr_counts_by_species(df::DataFrame; output_path::String="", max_species::Int=50)
    df_sorted = sort(df, [:clade, :n_acrr], rev=[false, true])

    if nrow(df_sorted) > max_species
        df_sorted = combine(groupby(df_sorted, :clade)) do sdf
            first(sort(sdf, :n_acrr, rev=true), min(5, nrow(sdf)))
        end
    end

    clades = unique(df_sorted.clade)
    clade_colors = Dict(c => Makie.wong_colors()[mod1(i, 7)] for (i, c) in enumerate(clades))

    fig = Figure(size=(1000, max(600, nrow(df_sorted) * 15)))
    ax = Axis(fig[1, 1],
        xlabel = "Number of ACRR pairs",
        ylabel = "Species",
        title = "Absolute Concentration Ratio Robustness by Species",
        yticks = (1:nrow(df_sorted), df_sorted.species_id),
        yticklabelsize = 8
    )

    colors = [clade_colors[c] for c in df_sorted.clade]
    barplot!(ax, 1:nrow(df_sorted), df_sorted.n_acrr,
        direction=:x, color=colors)

    if !isempty(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

"""
    plot_jaccard_vs_phylodist(jaccard_matrix, species_order) -> Figure

Scatter plot of Jaccard similarity vs pairwise index (proxy for phylo distance).
Note: For real phylogenetic distance, need to parse tree branch lengths.
"""
function plot_jaccard_scatter(jaccard_matrix::Matrix, df::DataFrame; output_path::String="")
    n = size(jaccard_matrix, 1)

    # Extract upper triangle (excluding diagonal)
    jaccard_values = Float64[]
    clade_same = Bool[]

    species_to_clade = Dict(row.species_id => row.clade for row in eachrow(df))
    species_order = df.species_id

    for i in 1:n
        for j in (i+1):n
            push!(jaccard_values, jaccard_matrix[i, j])
            push!(clade_same, get(species_to_clade, species_order[i], "") ==
                              get(species_to_clade, species_order[j], ""))
        end
    end

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1],
        xlabel = "Species pair index",
        ylabel = "Jaccard similarity (ACR sets)",
        title = "ACR Set Similarity Between Species Pairs"
    )

    # Color by same/different clade
    colors = [same ? :blue : :red for same in clade_same]
    scatter!(ax, 1:length(jaccard_values), jaccard_values,
        color=colors, markersize=3, alpha=0.5)

    # Add legend
    Legend(fig[1, 2],
        [MarkerElement(color=:blue, marker=:circle),
         MarkerElement(color=:red, marker=:circle)],
        ["Same clade", "Different clade"])

    if !isempty(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

"""
    plot_acr_heatmap(acr_matrix, species_order) -> Figure

Heatmap of ACR metabolite presence/absence across species.
"""
function plot_acr_heatmap(acr_matrix::DataFrame; output_path::String="", top_n_mets::Int=50)
    # Get conservation stats and filter to top N metabolites
    conservation = get_metabolite_conservation(acr_matrix)
    top_mets = first(conservation.metabolite, top_n_mets)

    # Filter matrix
    filtered = filter(row -> row.metabolite in top_mets, acr_matrix)

    # Convert to matrix for plotting
    mat = Matrix(filtered[:, 2:end])

    fig = Figure(size=(1400, 800))
    ax = Axis(fig[1, 1],
        xlabel = "Species",
        ylabel = "ACR Metabolite",
        title = "ACR Metabolite Presence/Absence (Top $top_n_mets)",
        xticklabelrotation = π/2,
        xticklabelsize = 6,
        yticklabelsize = 8
    )

    heatmap!(ax, mat', colormap=[:white, :darkblue])

    if !isempty(output_path)
        save(output_path, fig)
        println("Saved: $output_path")
    end

    return fig
end

# =============================================================================
# Clade Summary Statistics
# =============================================================================

"""
    calculate_clade_summary(df) -> DataFrame

Calculate summary statistics by clade.
"""
function calculate_clade_summary(df::DataFrame)
    summary = combine(groupby(df, :clade),
        :n_acr => mean => :mean_acr,
        :n_acr => std => :std_acr,
        :n_acr => median => :median_acr,
        :n_acrr => mean => :mean_acrr,
        :n_acrr => median => :median_acrr,
        :giant_proportion => mean => :mean_giant_prop,
        :giant_proportion => std => :std_giant_prop,
        :giant_proportion => median => :median_giant_prop,
        :largest_module_size => mean => :mean_giant_size,
        :largest_module_size => maximum => :max_giant_size,
        nrow => :n_species
    )
    sort!(summary, :mean_giant_prop, rev=true)
    return summary
end

# =============================================================================
# Main Analysis
# =============================================================================

function main()
    println("=" ^ 60)
    println("Kinetic Module & ACR Phylogenetic Analysis")
    println("=" ^ 60)

    # Create output directory
    mkpath(OUTPUT_DIR)
    println("Output directory: $OUTPUT_DIR")

    # 1. Load metadata
    println("\n1. Loading species metadata...")
    metadata = load_species_metadata(METADATA_FILE)
    println("   Loaded $(nrow(metadata)) species from metadata")

    # 2. Load kinetic results
    println("\n2. Loading kinetic analysis results...")
    kinetic_df, acr_dict = load_all_kinetic_results(RESULTS_DIR)

    # 3. Merge with metadata to get clade assignments
    println("\n3. Merging with clade assignments...")
    df = leftjoin(kinetic_df, metadata[:, [:species_id, :clade]], on=:species_id)

    # Fill missing clades
    df.clade = coalesce.(df.clade, "Unknown")

    println("   Species with clade assignments: $(count(!ismissing, df.clade))")
    println("   Unique clades: $(length(unique(df.clade)))")

    # 4. Save species summary
    println("\n4. Saving species kinetic summary...")
    summary_path = joinpath(OUTPUT_DIR, "species_kinetic_summary.csv")
    CSV.write(summary_path, df)
    println("   Saved: $summary_path")

    # 5. Calculate and save clade summary
    println("\n5. Calculating clade summary statistics...")
    clade_summary = calculate_clade_summary(df)
    clade_path = joinpath(OUTPUT_DIR, "clade_summary.csv")
    CSV.write(clade_path, clade_summary)
    println("   Saved: $clade_path")

    # Print clade summary
    println("\n   Clade Summary (sorted by mean giant module proportion):")
    println("   " * "-"^80)
    for row in eachrow(clade_summary)
        println("   $(rpad(row.clade, 35)) | n=$(lpad(row.n_species, 3)) | ACR=$(lpad(round(row.mean_acr, digits=1), 5)) | Giant=$(lpad(round(row.mean_giant_prop*100, digits=1), 5))%")
    end

    # 6. Build ACR presence matrix
    println("\n6. Building ACR presence/absence matrix...")
    species_order = df.species_id
    acr_matrix = build_acr_presence_matrix(acr_dict, species_order)
    acr_matrix_path = joinpath(OUTPUT_DIR, "acr_presence_matrix.csv")
    CSV.write(acr_matrix_path, acr_matrix)
    println("   Saved: $acr_matrix_path")

    # 7. Calculate metabolite conservation
    println("\n7. Calculating metabolite conservation...")
    conservation = get_metabolite_conservation(acr_matrix)

    # Add decoded (human-readable) IDs
    conservation.decoded_id = [decode_id(m) for m in conservation.metabolite]

    # Reorder columns
    select!(conservation, :metabolite, :decoded_id, :n_species_acr, :proportion)

    conservation_path = joinpath(OUTPUT_DIR, "acr_metabolite_conservation.csv")
    CSV.write(conservation_path, conservation)
    println("   Saved: $conservation_path")

    # Print top conserved metabolites
    println("\n   Top 20 most conserved ACR metabolites:")
    println("   " * "-"^80)
    for row in first(eachrow(conservation), 20)
        decoded = length(row.decoded_id) > 45 ? row.decoded_id[1:42] * "..." : row.decoded_id
        println("   $(rpad(decoded, 48)) | $(lpad(row.n_species_acr, 3)) species ($(lpad(round(row.proportion*100, digits=1), 5))%)")
    end

    # 8. Calculate Jaccard similarities
    println("\n8. Calculating pairwise Jaccard similarities...")
    jaccard_matrix = calculate_all_jaccard_similarities(acr_dict, species_order)

    # Save Jaccard matrix
    jaccard_df = DataFrame(jaccard_matrix, Symbol.(species_order))
    insertcols!(jaccard_df, 1, :species => species_order)
    jaccard_path = joinpath(OUTPUT_DIR, "jaccard_similarity_matrix.csv")
    CSV.write(jaccard_path, jaccard_df)
    println("   Saved: $jaccard_path")

    # 9. Generate plots
    println("\n9. Generating plots...")

    # Plot 1: Giant module by clade
    plot_giant_module_by_clade(df,
        output_path=joinpath(OUTPUT_DIR, "fig2_giant_module_by_clade.png"))

    # Plot 2: ACR counts
    plot_acr_counts_by_species(df,
        output_path=joinpath(OUTPUT_DIR, "fig3_acr_counts_by_species.png"))

    # Plot 3: ACRR counts
    plot_acrr_counts_by_species(df,
        output_path=joinpath(OUTPUT_DIR, "fig4_acrr_counts_by_species.png"))

    # Plot 4: Jaccard scatter
    plot_jaccard_scatter(jaccard_matrix, df,
        output_path=joinpath(OUTPUT_DIR, "fig6_jaccard_similarity.png"))

    # Plot 5: ACR heatmap
    plot_acr_heatmap(acr_matrix,
        output_path=joinpath(OUTPUT_DIR, "fig5_acr_heatmap.png"))

    println("\n" * "=" ^ 60)
    println("Analysis complete!")
    println("Results saved to: $OUTPUT_DIR")
    println("=" ^ 60)

    return df, acr_matrix, jaccard_matrix, clade_summary
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
