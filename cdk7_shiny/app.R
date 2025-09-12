# CDK7 Gigantic Heatmap Shiny App
library(shiny)
library(shinydashboard)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(tidyverse)
library(qs)
library(impute)
library(dplyr)
library(tibble)
library(grid)
library(seriation)
library(ggbeeswarm)

# Load data and extract parameter options
heatmap_data <- qread("../data/gigantic_heatmap_data.qs")

# Extract unique parameter values from data
gr_agents <- sort(unique(heatmap_data$baseline_deep_varstab_cor_df$agent))
gr_metrics <- sort(unique(heatmap_data$baseline_deep_varstab_cor_df$gr_metric))
crispr_screens <- sort(unique(heatmap_data$crispr_joined$screen))
gene_expr_conditions <- sort(unique(heatmap_data$deseq_res_acquired_data$treatment))
gene_expr_cell_lines <- sort(unique(heatmap_data$deseq_res_acquired_data$cell_line))
emm_contrasts <- sort(unique(heatmap_data$deseq_res_acquired_ykl$contrast))
emm_cell_lines <- sort(unique(heatmap_data$deseq_res_acquired_ykl$cell_line))

gene_id_gene_name_map <- qread("../data/ensembl_111_gene_id_name_map.qs")

# Load VST counts data for gene-specific visualization
counts_acquired_long <- read_csv("../deseq/cdk7_counts_acquired_long.csv.gz", show_col_types = FALSE)

# Available gene expression metrics
gene_expr_metrics <- c("signed_p", "log2FoldChange", "pvalue", "padj")

# Default selections (matching original analysis)
default_gr_agents <- c("YKL-5-124", "SY-5609")
default_gr_metrics <- c("GR50", "GRmax")
default_crispr_screens <- c("SYcrispri", "YKLcrispri", "SYcrispra", "YKLcrispra")
default_gene_expr_conditions <- c("Y1uM", "S100nM", "YKLr1", "YKLr2", "YKLr3", "SYr3", "YKLrTQr1", "YKLrTQr3")
default_gene_expr_cell_lines <- c("OVCAR4", "SNU8", "OVCAR8", "OVCAR3")
default_emm_contrasts <- c("YKLr - no")
default_emm_cell_lines <- c("all_cells")

# Recreate processing functions
prep_gigantic_hm <- function(
  gr_cor_agents = default_gr_agents,
  gr_cor_metrics = default_gr_metrics,
  crispr_screens = default_crispr_screens,
  gene_expression_conditions = default_gene_expr_conditions,
  gene_expression_metric = "signed_p",
  gene_expression_cell_lines = default_gene_expr_cell_lines,
  gene_expression_emm_contrasts = default_emm_contrasts,
  gene_expression_emm_cell_lines = default_emm_cell_lines,
  n_min = 2L
) {
  ge_metric_sym <- sym(gene_expression_metric)

  cor_df <- heatmap_data$baseline_deep_varstab_cor_df %>%
    filter(
      agent %in% gr_cor_agents,
      gr_metric %in% gr_cor_metrics
    )

  crispr_df <- heatmap_data$crispr_joined %>%
    filter(screen %in% crispr_screens)

  expression_df <- heatmap_data$deseq_res_acquired_data %>%
    rename(ensembl_gene_id = gene_id) %>%
    mutate(across(all_of(gene_expression_metric), \(x) replace_na(x, 0))) %>%
    filter(
      treatment %in% gene_expression_conditions,
      cell_line %in% gene_expression_cell_lines
    )

  emm_df <- heatmap_data$deseq_res_acquired_ykl %>%
    rename(ensembl_gene_id = gene_id) %>%
    mutate(across(all_of(gene_expression_metric), \(x) replace_na(x, 0))) %>%
    filter(
      contrast %in% gene_expression_emm_contrasts,
      cell_line %in% gene_expression_emm_cell_lines
    )

  genes_to_keep <- list(
    cor_df %>%
      filter(p < .01),
    crispr_df %>%
      filter(`Mann-Whitney p-value` < .01),
    expression_df %>%
      filter(padj < .01) %>%
      anti_join(
        heatmap_data$deseq_res_resistant_vs_treatment_classes %>%
          filter(combined_simple == "both_same"),
        by = c("ensembl_gene_id" = "gene_id", "cell_line", "treatment")
      ),
    emm_df %>%
      filter(padj < .01)
  ) %>%
    map("ensembl_gene_id") %>%
    {
      gene_counts <- table(unlist(.))
      names(gene_counts)[gene_counts >= n_min]
    }

  if (length(genes_to_keep) == 0) {
    return(NULL)
  }

  cor_df_formatted <- cor_df %>%
    select(agent, ensembl_gene_id, gr_metric, r) %>%
    pivot_wider(names_from = c(agent, gr_metric), values_from = r)

  crispr_df_formatted <- crispr_df %>%
    select(ensembl_gene_id, screen, signed_p) %>%
    pivot_wider(names_from = screen, values_from = signed_p)

  expression_df_formatted <- expression_df %>%
    select(ensembl_gene_id, cell_line, treatment, {{ge_metric_sym}}) %>%
    pivot_wider(names_from = c(treatment, cell_line), values_from = {{ge_metric_sym}})

  emm_df_formatted <- emm_df %>%
    select(ensembl_gene_id, cell_line, contrast, {{ge_metric_sym}}) %>%
    pivot_wider(names_from = c(contrast, cell_line), values_from = {{ge_metric_sym}})

  list(
    cor_df_formatted,
    crispr_df_formatted,
    expression_df_formatted,
    emm_df_formatted
  ) %>%
    map(
      \(x) filter(x, ensembl_gene_id %in% genes_to_keep)
    ) %>%
    reduce(full_join, by = "ensembl_gene_id")
}

scale_df <- function(df) {
  df %>%
    mutate(
      across(
        -ensembl_gene_id,
        \(x) scale(x, center = FALSE, scale = TRUE)[, 1]
      )
    )
}

cluster_fun_eucl_imp <- function(mat, features_in_rows = TRUE) {
  if (!features_in_rows) {
    mat <- t(mat)
  }
  mat_imp <- impute.knn(
    mat, rng.seed = 42, colmax = 90
  )[["data"]]
  if (!features_in_rows) {
    mat_imp <- t(mat_imp)
    mat <- t(mat)
  }
  dist_mat <- dist(mat_imp)
  clust <- hclust(dist_mat, method = "complete")
  reorder(clust, dist_mat, method = "OLO")
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "CDK7 Gigantic Heatmap"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Heatmap Parameters", tabName = "parameters", icon = icon("sliders")),

      br(),

      # GR Correlation Parameters
      h4("GR Correlation"),
      checkboxGroupInput(
        "gr_agents",
        "Agents:",
        choices = setNames(gr_agents, gr_agents),
        selected = default_gr_agents
      ),

      checkboxGroupInput(
        "gr_metrics",
        "Metrics:",
        choices = setNames(gr_metrics, gr_metrics),
        selected = default_gr_metrics
      ),

      br(),

      # CRISPR Parameters
      h4("CRISPR Screens"),
      checkboxGroupInput(
        "crispr_screens",
        "Screens:",
        choices = setNames(crispr_screens, crispr_screens),
        selected = default_crispr_screens
      ),

      br(),

      # Gene Expression Parameters
      h4("Gene Expression"),
      checkboxGroupInput(
        "gene_expr_conditions",
        "Conditions:",
        choices = setNames(gene_expr_conditions, gene_expr_conditions),
        selected = default_gene_expr_conditions
      ),

      selectInput(
        "gene_expr_metric",
        "Metric:",
        choices = setNames(gene_expr_metrics, gene_expr_metrics),
        selected = "signed_p"
      ),

      checkboxGroupInput(
        "gene_expr_cell_lines",
        "Cell Lines:",
        choices = setNames(gene_expr_cell_lines, gene_expr_cell_lines),
        selected = default_gene_expr_cell_lines
      ),

      br(),

      # EMM Parameters
      h4("EMM Contrasts"),
      checkboxGroupInput(
        "emm_contrasts",
        "Contrasts:",
        choices = setNames(emm_contrasts, emm_contrasts),
        selected = default_emm_contrasts
      ),

      checkboxGroupInput(
        "emm_cell_lines",
        "Cell Lines:",
        choices = setNames(emm_cell_lines, emm_cell_lines),
        selected = default_emm_cell_lines
      ),

      br(),

      # Filtering Parameters
      h4("Filtering"),
      numericInput(
        "n_min",
        "Minimum Gene Count:",
        value = 2,
        min = 1,
        max = 10
      ),

      br(),

      actionButton("generate", "Generate Heatmap", class = "btn-primary", width = "100%")
    )
  ),

  dashboardBody(
    fluidRow(
      box(
        status = "primary",
        width = 6,
        InteractiveComplexHeatmapOutput(
          "heatmap",
          title1 = "",
          height1 = "750px",
          compact = TRUE
        )
      ),
      box(
        status = "primary",
        width = 6,
        InteractiveComplexHeatmapOutput(
          "subheatmap",
          title1 = "",
          height1 = "750px",
          compact = TRUE
        )
      )
    ),
    fluidRow(
      box(
        title = "Gene VST Counts",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        uiOutput("gene_plot_ui")
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # Store selected gene ID and brushed selection
  selected_gene <- reactiveVal(NULL)
  brushed_row_indices <- reactiveVal(NULL)

  # Reactive heatmap data generation
  heatmap_matrix <- eventReactive(input$generate, {

    # Validate inputs
    validate(
      need(length(input$gr_agents) > 0, "Please select at least one GR agent"),
      need(length(input$gr_metrics) > 0, "Please select at least one GR metric"),
      need(length(input$crispr_screens) > 0, "Please select at least one CRISPR screen")
    )

    withProgress(message = "Generating heatmap...", {

      setProgress(0.1, detail = "Preparing data...")

      # Generate heatmap data
      hm_data <- prep_gigantic_hm(
        gr_cor_agents = input$gr_agents,
        gr_cor_metrics = input$gr_metrics,
        crispr_screens = input$crispr_screens,
        gene_expression_conditions = input$gene_expr_conditions,
        gene_expression_metric = input$gene_expr_metric,
        gene_expression_cell_lines = input$gene_expr_cell_lines,
        gene_expression_emm_contrasts = input$emm_contrasts,
        gene_expression_emm_cell_lines = input$emm_cell_lines,
        n_min = input$n_min
      )

      validate(
        need(!is.null(hm_data), "No genes pass the filtering criteria. Try lowering n_min or adjusting other parameters."),
        need(nrow(hm_data) > 0, "No genes found with current parameters.")
      )

      setProgress(0.2, detail = "Scaling data...")

      # Scale data and convert to matrix
      hm_mat <- hm_data %>%
        scale_df() %>%
        column_to_rownames("ensembl_gene_id") %>%
        as.matrix()

      validate(
        need(nrow(hm_mat) > 1, "Need at least 2 genes for clustering"),
        need(ncol(hm_mat) > 0, "No data columns available")
      )

      setProgress(0.4, detail = "Clustering...")
1
      # Perform clustering
      hm_row_clust <- cluster_fun_eucl_imp(hm_mat)

      setProgress(0.8, detail = "Creating heatmap...")

      # Create heatmap
      hm <- Heatmap(
        hm_mat,
        name = "Scaled\nsignal",
        na_col = "grey",
        cluster_rows = hm_row_clust,
        row_labels = gene_id_gene_name_map[rownames(hm_mat)],
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          title = "Scaled\nsignal",
          title_position = "topcenter",
          legend_height = unit(4, "cm")
        )
      )

      setProgress(1, detail = "Complete!")

      # Clear previous selections when new heatmap is generated
      brushed_row_indices(NULL)
      selected_gene(NULL)

      hm
    })
  })

  # Render main interactive heatmap with brush action
  observe({
    hm <- heatmap_matrix()
    if (!is.null(hm)) {
      makeInteractiveComplexHeatmap(
        input, output, session, hm, "heatmap",
        brush_action = function(df, output) {
          if (!is.null(df) && nrow(df) > 0) {
            # Get row indices from brush selection
            row_indices <- unique(unlist(df$row_index))
            brushed_row_indices(row_indices)
          } else {
            brushed_row_indices(NULL)
          }
        }
      )
    }
  })

  # Create reactive second heatmap from brushed selection
  subheatmap_matrix <- reactive({
    hm <- heatmap_matrix()
    row_indices <- brushed_row_indices()

    if (!is.null(hm) && !is.null(row_indices) && length(row_indices) > 0) {
      sub_mat <- hm@matrix[row_indices, , drop = FALSE]

      # Only create heatmap if we have enough genes
      if (nrow(sub_mat) > 1) {
        # Perform clustering on subset
        sub_row_clust <- cluster_fun_eucl_imp(sub_mat)

        # Create sub-heatmap
        Heatmap(
          sub_mat,
          name = "Scaled\nsignal",
          na_col = "grey",
          cluster_rows = sub_row_clust,
          row_labels = gene_id_gene_name_map[rownames(sub_mat)],
          cluster_columns = FALSE,
          show_row_names = TRUE,
          show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 10),
          heatmap_legend_param = list(
            title = "Scaled\nsignal",
            title_position = "topcenter",
            legend_height = unit(4, "cm")
          )
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  # Render second interactive heatmap with click action
  observe({
    sub_hm <- subheatmap_matrix()
    if (!is.null(sub_hm)) {
      makeInteractiveComplexHeatmap(
        input, output, session, sub_hm, "subheatmap",
        click_action = function(df, output) {
          if (!is.null(df) && nrow(df) > 0) {
            # Get the first selected gene ID from the sub-heatmap row
            gene_id <- rownames(sub_hm@matrix)[df[1, "row_index"]]
            selected_gene(gene_id)
          }
        }
      )
    }
  })

  # Initialize with default heatmap on startup
  observe({
    updateActionButton(session, "generate", label = "Generate Heatmap")
  })

  # Render gene-specific VST count plot
  output$gene_plot_ui <- renderUI({
    gene_id <- selected_gene()

    if (is.null(gene_id)) {
      div(
        style = "text-align: center; padding: 50px;",
        h4("Click on a gene in the heatmap to view VST counts"),
        p("Select a gene from the interactive heatmap above to see its expression pattern across samples.")
      )
    } else {
      # Filter data for selected gene
      gene_data <- counts_acquired_long %>%
        filter(
          type == "vst",
          ensembl_gene_id == gene_id
        )

      if (nrow(gene_data) == 0) {
        div(
          style = "text-align: center; padding: 50px;",
          h4("No VST count data available"),
          p(paste("No data found for gene:", gene_id))
        )
      } else {
        # Get gene name for title
        gene_name <- gene_id_gene_name_map[gene_id]
        if (is.na(gene_name)) gene_name <- gene_id

        # Create the plot
        p <- gene_data %>%
          ggplot(
            aes(
              x = count,
              y = interaction(cell_line, resistant_overall),
              color = drug_treatment,
              shape = resistant_batch
            )
          ) +
          geom_quasirandom() +
          labs(
            title = paste("VST Counts for", gene_name, paste0("(", gene_id, ")")),
            x = "VST Count",
            y = "Cell Line × Resistant Status",
            color = "Drug Treatment",
            shape = "Resistant Batch"
          ) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5)
          )

        renderPlot({
          p
        }, height = 400)
      }
    }
  })
}

# Run app
shinyApp(ui = ui, server = server)
