using JLD2
using CSV
using DataFrames
using Base.Threads
using Printf
using Statistics
using SparseArrays
# COCOA package is required for proper deserialization of ConcordanceResults
using COCOA

# ConcordanceType enum values from COCOA.data_structures:
# None = 0
# Concordant = 1
# Trivially_concordant = 2
# Balanced = 3
# Trivially_balanced = 4
# Not_present = 5 (our addition for the matrix)

# --- SETUP PATHS ---
# JLD2 no-split results: extract prpd_models/no_split.tar.gz from Zenodo,
# run kinetic analysis, then point WORK_DIR at the resulting JLD2 files.
const REPO_ROOT         = joinpath(@__DIR__, "..", "..")
const WORK_DIR          = get(ENV, "JLD2_NO_SPLIT_DIR", joinpath(REPO_ROOT, "data", "jld2", "no_split"))
const METADATA_FILE     = joinpath(REPO_ROOT, "data", "annotations", "cocoa_results_merged.csv")
const CONSERVATION_FILE = joinpath(REPO_ROOT, "results", "concordance", "core_complex_consistency.csv")
const OUTPUT_FILE       = joinpath(REPO_ROOT, "results", "concordance", "core_complex_property_matrix.csv")

println("\n", "="^80)
println("EXTRACTING ALL COMPLEX PROPERTY MATRIX")
println("="^80, "\n")

# --- IDENTIFY ALL UNIQUE COMPLEXES ---
println("Identifying all unique complexes across all species...")

# Match files to species first
files = readdir(WORK_DIR)
jld_files = filter(f -> endswith(f, ".jld2") && startswith(f, "kinetic_results_"), files)

# Load metadata to get species IDs
metadata_temp = CSV.read(METADATA_FILE, DataFrame)
all_species_ids = string.(metadata_temp.old_species_id)

# Sort species IDs by length (longest first) to avoid prefix matching issues
# e.g., Metschnikowia_matae_maris must be matched before Metschnikowia_matae
all_species_ids_sorted = sort(all_species_ids, by=length, rev=true)

# Create mapping: species_id -> filename
species_to_file_temp = Dict{String, String}()
for f in jld_files
    for id in all_species_ids_sorted
        if startswith(f, "kinetic_results_$(id)_")
            species_to_file_temp[id] = f
            break
        end
    end
end

# Collect all unique complexes
all_complexes_set = Set{String}()
println("  Scanning JLD2 files...")
for (i, (species_id, filename)) in enumerate(species_to_file_temp)
    if i % 50 == 0
        print(".")
        flush(stdout)
    end
    try
        filepath = joinpath(WORK_DIR, filename)
        data = JLD2.load(filepath)
        if haskey(data, "results")
            res = data["results"]
            for cid in res.complex_ids
                push!(all_complexes_set, string(cid))
            end
        end
    catch e
        # Skip files with errors
    end
end
println()

all_complexes = sort(collect(all_complexes_set))
n_complexes = length(all_complexes)
println("  Total unique complexes: ", n_complexes, "\n")

# --- LOAD SPECIES METADATA ---
println("Loading species metadata...")
metadata = CSV.read(METADATA_FILE, DataFrame)

# Create ordered list of all species (yeasts first, then outgroups)
all_species = DataFrame[]

# Yeasts (not outgroup)
yeasts = metadata[metadata."Major clade" .!= "Outgroup", :]
push!(all_species, yeasts)

# Outgroups
outgroups = metadata[metadata."Major clade" .== "Outgroup", :]
push!(all_species, outgroups)

species_list = vcat(all_species...)
species_ids = string.(species_list.old_species_id)
species_names = species_list."Species name"
n_species = length(species_ids)

println("  Total species: ", n_species)
println("    Yeasts: ", nrow(yeasts))
println("    Outgroups: ", nrow(outgroups), "\n")

# --- MATCH FILES TO SPECIES ---
println("Matching JLD2 files to species...")
files = readdir(WORK_DIR)
jld_files = filter(f -> endswith(f, ".jld2") && startswith(f, "kinetic_results_"), files)

# Sort species IDs by length (longest first) to avoid prefix matching issues
species_ids_sorted = sort(species_ids, by=length, rev=true)

# Create mapping: species_id -> filename
species_to_file = Dict{String, String}()

for f in jld_files
    for id in species_ids_sorted
        if startswith(f, "kinetic_results_$(id)_")
            species_to_file[id] = f
            break
        end
    end
end

println("  Matched files: ", length(species_to_file), " / ", n_species, "\n")

# --- INITIALIZE PROPERTY MATRIX ---
println("Initializing property matrix...")
# Matrix: rows = complexes, cols = species
# Values match ConcordanceType enum:
#   0 = None (singleton)
#   1 = Concordant
#   2 = Trivially_concordant
#   3 = Balanced
#   4 = Trivially_balanced
#   5 = Not present (our addition)
property_matrix = fill(5, n_complexes, n_species)  # Default to "Not present"
println("  Dimensions: ", n_complexes, " complexes × ", n_species, " species\n")

# --- EXTRACT PROPERTIES FROM JLD2 FILES ---
println("Extracting properties from JLD2 files...")
print("  Progress: ")
flush(stdout)

# Create mapping from complex ID to row index for fast lookup
complex_to_idx = Dict{String, Int}()
for (idx, c) in enumerate(all_complexes)
    complex_to_idx[string(c)] = idx
end

"""
Determine concordance type for a complex based on concordance_matrix.

For complexes in concordant modules (module_id > 0):
- If ANY relationship with module members is Concordant (1), return 1
- If ALL relationships are Trivially_concordant (2), return 2

For balanced complexes (module_id == 0):
- Check diagonal of concordance_matrix for Balanced (3) or Trivially_balanced (4)

For singletons (module_id == -1):
- Return None (0)
"""
function determine_concordance_type(res, complex_idx::Int)::Int
    module_id = res.concordance_modules[complex_idx]
    cm = res.concordance_matrix

    if module_id == -1
        # Singleton - None
        return 0
    elseif module_id == 0
        # Balanced - check diagonal for Balanced (3) vs Trivially_balanced (4)
        diag_val = cm[complex_idx, complex_idx]
        if diag_val == 4
            return 4  # Trivially_balanced
        else
            return 3  # Balanced (default for balanced module)
        end
    else
        # Concordant module - find concordance type from relationships
        # Find other members of the same module
        same_module_indices = findall(==(module_id), res.concordance_modules)

        has_concordant = false
        has_trivially_concordant = false

        for other_idx in same_module_indices
            if other_idx == complex_idx
                continue
            end
            # Get concordance type (ensure canonical order for upper triangular)
            i, j = min(complex_idx, other_idx), max(complex_idx, other_idx)
            conc_type = cm[i, j]

            if conc_type == 1
                has_concordant = true
            elseif conc_type == 2
                has_trivially_concordant = true
            end
        end

        # If any relationship is Concordant (1), the complex is Concordant
        # Otherwise, if all relationships are Trivially_concordant (2), use that
        if has_concordant
            return 1  # Concordant
        elseif has_trivially_concordant
            return 2  # Trivially_concordant
        else
            # Fallback (shouldn't happen for well-formed data)
            return 1  # Default to Concordant
        end
    end
end

n_processed = 0
for (species_idx, species_id) in enumerate(species_ids)
    if !haskey(species_to_file, species_id)
        continue  # No file for this species
    end

    filename = species_to_file[species_id]
    filepath = joinpath(WORK_DIR, filename)

    try
        # Load JLD2 file
        data = JLD2.load(filepath)
        if !haskey(data, "results")
            continue
        end
        res = data["results"]

        # Extract properties for each complex using full concordance type
        for (complex_idx, complex_id) in enumerate(res.complex_ids)
            complex_str = string(complex_id)

            # All complexes should be in our mapping
            if haskey(complex_to_idx, complex_str)
                row_idx = complex_to_idx[complex_str]

                # Determine concordance type from concordance_matrix
                property = determine_concordance_type(res, complex_idx)
                property_matrix[row_idx, species_idx] = property
            end
        end

        global n_processed += 1
        if n_processed % 50 == 0
            print(".")
            flush(stdout)
        end

    catch e
        println("\n  Error reading ", filename, ": ", e)
    end
end

println("\n  Processed files: ", n_processed, "\n")

# --- CHECK COMPLETENESS ---
println("Checking matrix completeness...")
n_not_present = sum(property_matrix .== 5)  # 5 = Not present
pct_present = 100 * (1 - n_not_present / (n_complexes * n_species))
println("  Not present: ", n_not_present, " / ", n_complexes * n_species,
        " (", round(pct_present, digits=2), "% present)\n")

# --- SAVE MATRIX ---
println("Saving property matrix...")

# Create DataFrame with proper column names
# Column 1: Complex ID
# Columns 2+: Species (using species_id)
col_names = ["Complex"; species_ids]
output_df = DataFrame(property_matrix, Symbol.(species_ids))
insertcols!(output_df, 1, :Complex => all_complexes)

# Save to CSV
mkpath(dirname(OUTPUT_FILE))
CSV.write(OUTPUT_FILE, output_df)

println("  Saved to: ", OUTPUT_FILE, "\n")

# --- SUMMARY STATISTICS ---
println("="^80)
println("SUMMARY")
println("="^80, "\n")

println("Matrix dimensions:")
println("  Complexes: ", n_complexes)
println("  Species: ", n_species)
println("  Total cells: ", n_complexes * n_species, "\n")

println("Property distribution (ConcordanceType enum values):")
n_none = sum(property_matrix .== 0)
n_concordant = sum(property_matrix .== 1)
n_trivially_concordant = sum(property_matrix .== 2)
n_balanced = sum(property_matrix .== 3)
n_trivially_balanced = sum(property_matrix .== 4)
n_not_present = sum(property_matrix .== 5)
n_total = n_complexes * n_species

println("  0 - None (singleton): ", n_none, " (", round(100 * n_none / n_total, digits=2), "%)")
println("  1 - Concordant: ", n_concordant, " (", round(100 * n_concordant / n_total, digits=2), "%)")
println("  2 - Trivially_concordant: ", n_trivially_concordant, " (", round(100 * n_trivially_concordant / n_total, digits=2), "%)")
println("  3 - Balanced: ", n_balanced, " (", round(100 * n_balanced / n_total, digits=2), "%)")
println("  4 - Trivially_balanced: ", n_trivially_balanced, " (", round(100 * n_trivially_balanced / n_total, digits=2), "%)")
println("  5 - Not present: ", n_not_present, " (", round(100 * n_not_present / n_total, digits=2), "%)")
println("  Total cells: ", n_total, "\n")

# Summary by category
n_concordant_total = n_concordant + n_trivially_concordant
n_balanced_total = n_balanced + n_trivially_balanced
println("Category summary:")
println("  Concordant (total): ", n_concordant_total, " (", round(100 * n_concordant_total / n_total, digits=2), "%)")
println("  Balanced (total): ", n_balanced_total, " (", round(100 * n_balanced_total / n_total, digits=2), "%)")
println("  Singleton (None): ", n_none, " (", round(100 * n_none / n_total, digits=2), "%)")
println("  Not present: ", n_not_present, " (", round(100 * n_not_present / n_total, digits=2), "%)\n")

println("Presence by complex:")
presence_per_complex = [sum(property_matrix[i, :] .!= 5) for i in 1:n_complexes]  # 5 = Not present
println("  Mean species per complex: ", round(mean(presence_per_complex), digits=1))
println("  Min: ", minimum(presence_per_complex))
println("  Max: ", maximum(presence_per_complex), " (out of ", n_species, " species)\n")

println("="^80)
println("EXTRACTION COMPLETE")
println("="^80, "\n")
