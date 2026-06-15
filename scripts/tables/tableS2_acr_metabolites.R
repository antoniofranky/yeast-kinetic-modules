################################################################################
## CLASSIFY ACR METABOLITES BY TYPE USING PATHWAY ANNOTATIONS
################################################################################

library(tidyverse)
library(xml2)

library(here)
REPO_ROOT   <- here::here()
DATA_DIR    <- file.path(REPO_ROOT, "data")
RESULTS_DIR <- file.path(REPO_ROOT, "results")
# NOTE: Saccharomyces_cerevisiae.xml is from the Yeast-Species-GEMs dataset (Shen et al. 2018).
# Download from https://github.com/SysBioChalmers/Yeast-Species-GEMs and set the path below.
YEAST_GEM_XML <- file.path(DATA_DIR, "Saccharomyces_cerevisiae.xml")

cat("Loading Yeast-GEM model to extract pathway annotations...\n")

# Read SBML XML
model_xml <- read_xml(YEAST_GEM_XML)

# Extract reaction-to-pathway mappings from groups
groups <- xml_find_all(model_xml, ".//groups:group", xml_ns(model_xml))

reaction_pathways <- tibble()

for (group in groups) {
  pathway_name <- xml_attr(group, "name")

  # Find reactions in this group
  members <- xml_find_all(group, ".//groups:member", xml_ns(model_xml))

  for (member in members) {
    rxn_id <- xml_attr(member, "idRef")
    reaction_pathways <- bind_rows(reaction_pathways,
                                   tibble(reaction_id = rxn_id, pathway = pathway_name))
  }
}

cat("Extracted", nrow(reaction_pathways), "reaction-pathway associations\n")
cat("Unique pathways:", length(unique(reaction_pathways$pathway)), "\n\n")

# Load template with names
acr <- read_csv(file.path(RESULTS_DIR, "acr/metabolite_classification_template_with_names.csv"),
                show_col_types = FALSE)

cat("Classifying", nrow(acr), "ACR metabolites using pathway annotations...\n\n")

# Extract metabolite ID-to-name mapping from SBML
cat("Extracting metabolite names from SBML...\n")
species_xml <- xml_find_all(model_xml, ".//d1:species", xml_ns(model_xml))

sbml_met_names <- tibble()
for (species in species_xml) {
  sbml_met_names <- bind_rows(sbml_met_names,
                              tibble(
                                sbml_id = xml_attr(species, "id"),
                                sbml_name = xml_attr(species, "name")
                              ))
}

cat("Extracted", nrow(sbml_met_names), "metabolite names from SBML\n")

# Extract metabolite-reaction associations from SBML
cat("Extracting metabolite-reaction associations...\n")
reactions_xml <- xml_find_all(model_xml, ".//d1:reaction", xml_ns(model_xml))

met_reaction_map <- tibble()

for (rxn_node in reactions_xml) {
  rxn_id <- xml_attr(rxn_node, "id")

  # Get reactants and products
  species_refs <- xml_find_all(rxn_node, ".//d1:speciesReference[@species]", xml_ns(model_xml))
  for (spec in species_refs) {
    met_id <- xml_attr(spec, "species")
    met_reaction_map <- bind_rows(met_reaction_map,
                                  tibble(sbml_id = met_id, reaction_id = rxn_id))
  }
}

met_reaction_map <- met_reaction_map %>% distinct()

cat("Extracted", nrow(met_reaction_map), "metabolite-reaction associations\n\n")

# Join metabolites -> reactions -> pathways
cat("Mapping ACR metabolites to pathways...\n")

# Match ACR metabolites to SBML by name
# ACR names are like "ATP [cytoplasm]", SBML names should match exactly
acr_name_map <- acr %>%
  left_join(sbml_met_names, by = c("name" = "sbml_name")) %>%
  filter(!is.na(sbml_id))

cat("  Matched", nrow(acr_name_map), "of", nrow(acr), "ACR metabolites to SBML by exact name\n")

# Now join to reactions and pathways
met_pathways <- acr_name_map %>%
  left_join(met_reaction_map, by = "sbml_id", relationship = "many-to-many") %>%
  left_join(reaction_pathways, by = "reaction_id", relationship = "many-to-many")

# Ensure pathway column exists (in case no matches were found)
if (!"pathway" %in% colnames(met_pathways)) {
  met_pathways <- met_pathways %>%
    mutate(pathway = NA_character_)
}

# Filter to keep only rows with pathway annotations
met_pathways <- met_pathways %>%
  filter(!is.na(pathway))

cat("  Mapped", length(unique(met_pathways$metabolite_id)), "ACR metabolites to pathways\n\n")

# Keep pathway annotations as-is (don't combine)
acr_pathway_annotations <- met_pathways %>%
  select(metabolite_id, name, pathway, reaction_id) %>%
  arrange(metabolite_id, pathway)

cat("Mapped", length(unique(acr_pathway_annotations$metabolite_id)), "ACR metabolites to pathways\n\n")

# Save pathway annotations
write_csv(acr_pathway_annotations,
          file.path(RESULTS_DIR, "acr/acr_pathway_annotations.csv"))

cat("Saved pathway annotations to: ../results/acr_evolution/acr_pathway_annotations.csv\n\n")

# For classification: assign primary pathway (most common pathway for each metabolite)
primary_pathways <- met_pathways %>%
  group_by(metabolite_id, pathway) %>%
  summarise(n_reactions = n(), .groups = "drop") %>%
  group_by(metabolite_id) %>%
  slice_max(n_reactions, n = 1, with_ties = FALSE) %>%
  select(metabolite_id, primary_pathway = pathway)

# Add primary pathway to ACR data (keep original pathway names, no broad categories)
acr_classified <- acr %>%
  left_join(primary_pathways, by = "metabolite_id")

# Summary by pathway
cat("Classification summary by pathway:\n")
pathway_summary <- acr_classified %>%
  filter(!is.na(primary_pathway)) %>%
  group_by(primary_pathway, conservation_class) %>%
  summarise(
    n = n(),
    mean_prevalence = mean(prevalence),
    .groups = "drop"
  ) %>%
  arrange(primary_pathway, desc(conservation_class))

print(pathway_summary, n = 100)

cat("\n\nOverall by pathway:\n")
overall_pathway <- acr_classified %>%
  filter(!is.na(primary_pathway)) %>%
  group_by(primary_pathway) %>%
  summarise(
    n = n(),
    pct = 100 * n() / sum(!is.na(acr_classified$primary_pathway)),
    mean_prevalence = mean(prevalence),
    n_core = sum(conservation_class == "Core (≥90%)"),
    pct_core = 100 * sum(conservation_class == "Core (≥90%)") / n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_prevalence))

print(overall_pathway, width = 120)

# Save classified data with original pathways
write_csv(acr_classified,
          file.path(RESULTS_DIR, "acr/metabolite_classification_complete.csv"))

cat("\n\nSaved classified metabolites to:")
cat("\n  ../results/acr_evolution/metabolite_classification_complete.csv\n")

# Also create a clean version for thesis
acr_thesis <- acr_classified %>%
  select(metabolite_id, name, compartment, formula, prevalence, conservation_class,
         primary_pathway)

write_csv(acr_thesis,
          file.path(RESULTS_DIR, "acr/acr_metabolites_thesis_table.csv"))

cat("  ../results/acr_evolution/acr_metabolites_thesis_table.csv\n\n")

cat("Done!\n")

