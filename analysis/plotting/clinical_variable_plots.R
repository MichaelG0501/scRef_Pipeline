####################
# Moved from: analysis/plot_clinical.R
# Reorganized as part of analysis/ restructuring
####################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr) 

library(dplyr)
library(readxl)

data <- read_excel("/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Concise_Summary_EAC_Ref.xlsx",
                   sheet = 3, skip = 1) %>%
  mutate(orig.ident = paste(Author, Year, `Sample Name`, sep = "_"))

cols_to_remove <- intersect(colnames(data), colnames(tmdata_all@meta.data))
cols_to_remove <- setdiff(cols_to_remove, "orig.ident")
tmdata_all@meta.data[, cols_to_remove] <- NULL

tmdata_all@meta.data <- tmdata_all@meta.data %>%
  left_join(data, by = "orig.ident")
cell_names <- colnames(tmdata_all)
rownames(tmdata_all@meta.data) <- cell_names


plot_all <- function(studyname = NULL) {
  
  # --- 1. Data Preparation ---
  df <- tmdata_all@meta.data
  
  # Check if we are filtering for a specific study
  if (!is.null(studyname)) {
    df <- df %>% filter(study == studyname)
    show_study_details <- FALSE # Flag: Hide detailed study list
    file_suffix <- studyname
  } else {
    studyname <- "All"
    show_study_details <- TRUE  # Flag: Show detailed study list
    file_suffix <- "All"
  }
  
  # Ensure numeric and factors are clean
  df$Age <- as.numeric(df$Age) 
  df$Age_Group <- ifelse(df$Age > 60, ">60", "<=60")
  
  # Standardize Gender
  df$Gender <- ifelse(df$Gender %in% c("F", "female", "Female"), "Female", "Male")
  
  # Factorize Treatment and Response to control order
  df$Treatment <- factor(df$Treatment, levels = c("Tx-naïve", "Post"))
  df$`Clinical response` <- factor(ifelse(df$`Clinical response` == "R", "Responder", "Nonresponder"), 
                                   levels = c("Responder", "Nonresponder"))
  
  # Define Colors (Ensure manual_names is available in your environment)
  manual_cols <- setNames(
    c("#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "grey80"),
    manual_names 
  )
  
  # --- 2. Enhanced Plotting Function ---
  plot_prop <- function(data, group_var, cell_state_var = "manual_state", 
                        study_col = "study", 
                        total_studies_N = 8, 
                        filter_expr = NULL, title = NULL) {
    
    # A. Filter Data
    if (!is.null(filter_expr)) {
      data <- data %>% filter(!!rlang::parse_expr(filter_expr))
    }
    data <- data %>% filter(!is.na(.data[[group_var]]))
    
    # B. Calculate Stats & Generate Labels
    stats_data <- data %>%
      group_by(.data[[group_var]]) %>%
      summarise(
        n_cells = n(),
        n_studies_present = n_distinct(.data[[study_col]]),
        studies_list = paste(sort(unique(.data[[study_col]])), collapse = "\n"),
        .groups = 'drop'
      ) %>%
      mutate(
        # CONDITIONAL LABEL LOGIC
        label_text = if (show_study_details) {
          # Complex Label: N + (x/8) + List of names
          header <- paste0("N=", comma(n_cells), "\n(", n_studies_present, "/", total_studies_N, ")")
          paste0(header, "\n\n", studies_list)
        } else {
          # Simple Label: N only (for single study plots)
          paste0("N=", comma(n_cells))
        }
      ) %>%
      mutate(
        # Count lines to estimate vertical space needed
        line_count = str_count(label_text, "\n") + 1
      )
    
    # C. Calculate Cell State Proportions
    plot_data <- data %>%
      group_by(.data[[group_var]], .data[[cell_state_var]]) %>%
      summarise(n = n(), .groups = 'drop') %>%
      group_by(.data[[group_var]]) %>%
      mutate(freq = n / sum(n))
    
    # D. Determine Dynamic Y-Axis Expansion
    max_lines <- max(stats_data$line_count)
    # If showing details, expand more. If simple, use small fixed expansion.
    expansion_mult <- if (show_study_details) {
      0.04 + (max_lines * 0.04) 
    } else {
      0.15 
    }
    
    # E. Plotting
    ggplot(plot_data, aes(x = .data[[group_var]], y = freq, fill = .data[[cell_state_var]])) +
      geom_bar(stat = "identity", position = "fill", width = 0.7) +
      
      # Text Annotation
      geom_text(data = stats_data, 
                aes(x = .data[[group_var]], y = 1.02, label = label_text), 
                inherit.aes = FALSE, 
                vjust = 0,      
                size = 3,       
                lineheight = 0.9, 
                color = "black") +
      
      # Formatting
      scale_y_continuous(labels = scales::percent, 
                         expand = expansion(mult = c(0, expansion_mult))) + 
      scale_fill_manual(values = manual_cols, drop = FALSE) +
      theme_minimal(base_size = 16) + 
      labs(x = NULL, y = "Proportion", title = title %||% group_var, fill = "Cell State") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", 
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 18)
      )
  }
  
  # --- 3. Generate Individual Plots ---
  p1 <- plot_prop(df, "Gender", title = "Gender")
  p2 <- plot_prop(df, "Age_Group", title = "Age (>60)")
  p3 <- plot_prop(df, "Race/Ethnicity", title = "Race")
  
  p4 <- plot_prop(df, "Tumor Location", title = "Location")
  p5 <- plot_prop(df, "Tumor Type", title = "Type")
  p6 <- plot_prop(df, "Grade (Differentiation)", title = "Grade")
  
  p7 <- plot_prop(df, "Stage AJCC", title = "Stage (All)")
  p8 <- plot_prop(df, "Stage AJCC", filter_expr = "Treatment == 'Tx-naïve'", title = "Stage (Tx-Naive)")
  p9 <- plot_prop(df, "Stage AJCC", filter_expr = "Treatment == 'Post'", title = "Stage (Post-Tx)")
  
  p10 <- plot_prop(df, "Technology", title = "Technology")
  p11 <- plot_prop(df, "Treatment", title = "Tx Timing")
  
  p12 <- plot_prop(df, "Clinical response", filter_expr = "Treatment == 'Tx-naïve'", title = "Response (Predictive)")
  p13 <- plot_prop(df, "Clinical response", filter_expr = "Treatment == 'Post'", title = "Response (Post-Tx)")
  
  
  # --- 4. Generate PDF Report ---
  pdf_filename <- paste0("CellState_clinical_", file_suffix, ".pdf")
  
  pdf(pdf_filename, width = 14, height = 16)
  
  # Page 1: Demographics
  print(
    ((p1 + p2) / p3) + 
      plot_layout(heights = c(1, 1.2)) + 
      plot_annotation(title = paste("Demographics Overview -", studyname))
  )
  
  # Page 2: Tumor Characteristics
  print(
    ((p4 + p5) / p6) + 
      plot_layout(heights = c(1, 1.2)) + 
      plot_annotation(title = paste("Tumor Characteristics -", studyname))
  )
  
  # Page 3: Staging
  print(
    (p7 / p8 / p9) + 
      plot_annotation(title = paste("Staging Analysis -", studyname))
  )
  
  # Page 4: Technology & Treatment
  print(
    (p10 / p11) + 
      plot_layout(heights = c(1.5, 1)) + 
      plot_annotation(title = paste("Technology & Treatment -", studyname))
  )
  
  # Page 5: Response
  print(
    (p12 + p13) + 
      plot_annotation(title = paste("Clinical Response -", studyname))
  )
  
  dev.off()
  
  message(paste("PDF generated:", pdf_filename))
}

plot_all()


###################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(stringr)

# --- 1. Data Preparation ---
df <- tmdata_all@meta.data

# Clean Numeric/Factor Variables
df$Age <- as.numeric(df$Age) 
df$Age_Group <- ifelse(df$Age > 60, ">60", "<=60")

# Standardize Gender
df$Gender <- ifelse(df$Gender %in% c("F", "female", "Female"), "Female", "Male")

# Factorize Treatment and Response
df$Treatment <- factor(df$Treatment, levels = c("Tx-naïve", "Post"))
df$`Clinical response` <- factor(ifelse(df$`Clinical response` == "R", "Responder", "Nonresponder"), 
                                 levels = c("Responder", "Nonresponder"))

# Define Colors
manual_cols <- setNames(
  c("#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "grey80"),
  manual_names 
)

# --- 2. Faceted Plotting Function ---
plot_variable_facet <- function(data, group_var, cell_state_var = "manual_state", 
                                study_col = "study", 
                                filter_expr = NULL, title = NULL) {
  
  # A. Filter Data (Global Filter for this plot)
  if (!is.null(filter_expr)) {
    data <- data %>% filter(!!rlang::parse_expr(filter_expr))
  }
  
  # B. Filter NAs for the specific variable
  # This ensures studies with NO data for this variable disappear from the facet
  data <- data %>% filter(!is.na(.data[[group_var]]))
  
  # If no data remains after filtering, return NULL (skip page)
  if (nrow(data) == 0) return(NULL)
  
  # C. Create Facet Labels with "N=" included
  # This calculates the total N for that study *filtered by the criteria*
  # e.g., "Alcindor_2025\n(N=450)"
  study_stats <- data %>%
    group_by(.data[[study_col]]) %>%
    mutate(study_label = paste0(.data[[study_col]], "\n(N=", comma(n()), ")")) %>%
    ungroup()
  
  # D. Calculate Proportions
  plot_data <- study_stats %>%
    group_by(study_label, .data[[group_var]], .data[[cell_state_var]]) %>%
    summarise(n = n(), .groups = 'drop') %>%
    group_by(study_label, .data[[group_var]]) %>%
    mutate(freq = n / sum(n))
  bar_totals <- plot_data %>%
    group_by(study_label, .data[[group_var]]) %>%
    summarise(total_n = sum(n), .groups = "drop")
  
  # E. Plot
  p <- ggplot(plot_data, aes(x = .data[[group_var]], y = freq, fill = .data[[cell_state_var]])) +
    geom_bar(stat = "identity", position = "fill", width = 0.8) +
    geom_text(data = bar_totals,
              aes(x = .data[[group_var]], y = 1.02, label = total_n),
              inherit.aes = FALSE,
              size = 3) +
    # Facet by Study (using the label with N)
    facet_wrap(~ study_label, scales = "free_x", ncol = 4) + 
    
    scale_y_continuous(labels = scales::percent) + 
    scale_fill_manual(values = manual_cols, drop = FALSE) +
    
    theme_bw(base_size = 14) + 
    labs(x = NULL, y = "Proportion", title = title %||% group_var, fill = "Cell State") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom", 
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95"), # Light grey header for study names
      strip.text = element_text(face = "bold", size = 11)
    )
  
  return(p)
}


# --- 3. Define the List of Plots ---
# We create a list of list to iterate through easily
plot_configs <- list(
  list(var = "Gender", title = "Gender Distribution by Study", filter = NULL),
  list(var = "Age_Group", title = "Age (>60) Distribution by Study", filter = NULL),
  list(var = "Race/Ethnicity", title = "Race Distribution by Study", filter = NULL),
  
  list(var = "Tumor Location", title = "Tumor Location by Study", filter = NULL),
  list(var = "Tumor Type", title = "Tumor Type by Study", filter = NULL),
  list(var = "Grade (Differentiation)", title = "Grade by Study", filter = NULL),
  
  list(var = "Stage AJCC", title = "Stage (All Samples) by Study", filter = NULL),
  list(var = "Stage AJCC", title = "Stage (Tx-Naive Only) by Study", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Stage AJCC", title = "Stage (Post-Tx Only) by Study", filter = "Treatment == 'Post'"),
  
  list(var = "Technology", title = "Technology Used by Study", filter = NULL),
  list(var = "Treatment", title = "Treatment Timing by Study", filter = NULL),
  
  list(var = "Clinical response", title = "Response (Predictive / Tx-Naive) by Study", filter = "Treatment == 'Tx-naïve'"),
  list(var = "Clinical response", title = "Response (Assessment / Post-Tx) by Study", filter = "Treatment == 'Post'")
)

# --- 4. Loop and Generate PDF ---

pdf("CellState_PerVariable_ByStudy.pdf", width = 16, height = 10) # Landscape for grid view

for (cfg in plot_configs) {
  
  message(paste("Plotting:", cfg$title))
  
  p <- plot_variable_facet(
    data = df, 
    group_var = cfg$var, 
    filter_expr = cfg$filter, 
    title = cfg$title
  )
  
  # Only print if the plot was generated (not NULL)
  if (!is.null(p)) {
    print(p)
  }
}

dev.off()
message("Done! Saved as CellState_PerVariable_ByStudy.pdf")

