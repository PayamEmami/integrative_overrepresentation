# pathway integration

This small script can be used to perform pathway-level data integration for a number of modalities where KEGG pathway exists. 

# Integrative Over-Representation Analysis with KEGG Pathways

## Description

The `do_pathway_overrepresentation` function is designed for performing integrative over-representation analysis using a set of pre-defined KEGG pathways. This powerful analysis tool allows you to investigate the enrichment of specific pathways in different data modalities, such as genes, metabolites, or other features. Whether you're exploring gene expression data, metabolomics data, or any other omics data, this function helps you identify pathways that are statistically enriched in your dataset.

## Required packages
You need to install the following libraries:

`KEGGREST`
`future.apply`
`foreach`
`doFuture`

## Usage

To use this function, follow these steps:

1. **Prepare your data**: Your data should be organized in a specific format with three essential columns: `pathway`, `hit`, and `name`. The `pathway` column should contain KEGG map IDs (e.g., map00010), with multiple pathways separated by commas or spaces if applicable. The `hit` column should be of logical type, where `TRUE` indicates that a feature is selected, and `FALSE` indicates it's not selected (e.g., based on a statistical test). The `name` column should contain unique names or IDs for each feature.

2. **Create a `references` list**: Create a named list of data frames, with each element representing a different data modality (e.g., genes, metabolites). Make sure that each data frame contains all quantified data, including both significant and non-significant values.

3. **Define `omics_type`**: Create a vector of character strings, with each element corresponding to the type of omics data in your `references` list (e.g., "ORTHOLOGY" for genes, "COMPOUND" for metabolites).

4. **Build a database**: Use the `create_database` function to create a special database, which you will pass as an argument to `do_pathway_overrepresentation`. See below

5. **Run the analysis**: Call the `do_pathway_overrepresentation` function with the appropriate arguments, including `references`, `omics_type`, and `database`.

## Parameters

- `references`: A named list of data frames, where each element represents a different data modality.
- `omics_type`: A vector of character strings indicating the type of omics data for each element in `references`.
- `database`: A special database created using the `create_database` function.
- `test_type`: Specify the statistical test to perform (Fisher or hypergeometric). Default is hypergeometric.
- `background_type`: Choose between local or global KEGG background for the analysis. Default is local.
- `weight_background_type`: Determine whether to use local or global KEGG data for weighting different omics. Default is global.
- `merged_analysis`: If `TRUE`, feature merging is performed instead of pathway weighting.
- `weight_type`: Select one of "unweighted," "overall," or "pathway" for weighting pathways differently.
- `pattern`: Specify the pattern for extracting KEGG IDs (leave as is unless you have specific requirements).
- `pathway_column`: Define the column in the data frames that contains pathway IDs.
- `hit_column`: Identify the column in the data frames that indicates hits.
- `name_column`: Specify the column in the data frames that contains feature names.

## Output

The function returns a list with pathway information for each omics data modality and the overlapping pathways.

## Details

- Each data frame in `references` should contain three columns: `pathway`, `hit`, and `name`.
- The `pathway` column should contain KEGG map IDs, with multiple pathways separated by commas or spaces.
- The `hit` column should be of logical type (TRUE or FALSE).
- The `references` parameter is a list of data frames, with each element representing a data modality.
- `omics_type` is a vector indicating the type of omics data for each element in `references` (e.g., "ORTHOLOGY" or "COMPOUND").
- `background_type` determines whether to use local or global KEGG background data.
- `weight_background_type` controls the weighting of omics data based on local or global KEGG data.
- `weight_type` specifies the weighting scheme for pathways.
- For the best results, consider using `weight_background_type` as "global."

## Example

```R
# Example data frames
metabolomics_data <- data.frame(
  name = c("metabo1", "metabo2", "metabo3"),
  pathway = c("map00930,map01100,map01120,map01220", "map01502", "map01502,map02010,map04978,map0523"),
  hit = c(FALSE, FALSE, TRUE)
)

transcript_data <- data.frame(
  name = c("gene1", "gene2", "gene3", "gene4"),
  pathway = c("map00930,map01100", "map01502,map00930", "map01502,map02010,map04978", "map01502,map02010"),
  hit = c(FALSE, FALSE, TRUE, TRUE)
)

mirna_data <- data.frame(
  name = c("m1", "m2", "m3", "m4"),
  pathway = c("map00930", "map01502,map00930", "map01502,map02010", "map01502,map02010"),
  hit = c(FALSE, FALSE, TRUE, TRUE)
)

# Create references list
references <- list(
  transcripts = transcript_data,
  metabolites = metabolomics_data,
  mirna = mirna_data
)

# Define omics_type
omics_type <- c("ORTHOLOGY", "COMPOUND", "ORTHOLOGY")

# Create a database
database <- create_database(references)

# Perform pathway over-representation analysis
results <- do_pathway_overrepresentation(
  references = references,
  omics_type = omics_type,
  database_input = database
)
```

# Extract KEGG Database Entries using create_database

## Description

The `extract_kegg_database_entries` function is designed to extract KEGG database entries, primarily used to keep track of the number of genes and compounds in the KEGG data.

## Usage

To use this function, follow these steps:

1. **Prepare your data**: Your data should be organized in a specific format with columns that include pathway IDs, hits, and compound or gene names.

2. **Create a `references` list**: Build a named list of data frames, where each element represents a different data modality (e.g., genes, metabolites). Ensure that each data frame contains the relevant information, including pathway IDs, hits, and names.

3. **Specify columns**: Define the column names in your data frames that contain pathway IDs (`pathway_column`), hits (`hit_column`), and compound or gene names (`name_column`).

4. **Set the number of cores**: Specify the number of CPU cores to use for the extraction process using the `ncores` parameter. Be cautious not to use an excessively high number of cores, as this may lead to a block in the connection to the KEGG database.

## Parameters

- `references`: A named list of data frames, where each element represents a different data modality.
- `pathway_column`: Identify the column in the data frames that contains pathway IDs.
- `hit_column`: Identify the column in the data frames that indicates hits.
- `name_column`: Specify the column in the data frames that contains feature names.
- `ncores`: Number of CPU cores to use for the extraction process.

## Output

The function returns a list of pathways with the associated entries information from the KEGG database.

## Details

- Be cautious with the `ncores` parameter, as setting it too high may lead to connection issues with the KEGG database. A reasonable number of cores, typically between 1 and 5, should suffice. If you encounter numerous warnings, consider reducing the number of cores.
- Keep in mind that some pathway maps may not match perfectly, so it's important to verify that any warnings do not result from missing pathways.

## Example

```R
metabolomics_data <- data.frame(name=c("metabo1","metabo2","metabo3"),
                                pathway=c("map00930,map01100,map01120,map01220","map01502","map01502,map02010,map04978,map0523"),
                                hit=c(FALSE,FALSE,TRUE))

transcript_data <- data.frame(name=c("gene1","gene2","gene3","gene4"),
                              pathway=c("map00930,map01100","map01502,map00930","map01502,map02010,map04978","map01502,map02010"),
                              hit=c(FALSE,FALSE,TRUE,TRUE))

mirna_data <- data.frame(name=c("m1","m2","m3","m4"),
                         pathway=c("map00930","map01502,map00930","map01502,map02010","map01502,map02010"),
                         hit=c(FALSE,FALSE,TRUE,TRUE))
references <- list(transcripts=transcript_data,metabolites=metabolomics_data,mirna=mirna_data)
database <- create_database(references)
```
