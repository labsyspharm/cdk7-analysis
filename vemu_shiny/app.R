library(shiny)
library(plotly)
library(tidyverse)
library(here)
library(qs)
library(ggbeeswarm)
library(data.table)
library(InteractiveComplexHeatmap)

theme_set(theme_minimal())

PATH_PREFIX = "vemu"

cluster_mat <- function(mat) {
  d <- dist(mat)
  set.seed(42)
  hclust(d, method = "average") %>%
    seriation:::reorder.hclust(dist = d, method = "olo")
}

find_scale_range <- function(x, percentile = 0.05) {
  abs_max <- quantile(x, c(percentile, 1 - percentile), names = FALSE, na.rm = TRUE) %>%
    abs() %>%
    max()
  c(-abs_max, abs_max)
}

cluster_df <- function(df, x, y, z) {
  # browser()
  mat <- df %>%
    dplyr::select({{x}}, {{y}}, {{z}}) %>%
    pivot_wider(
      names_from = {{x}},
      values_from = {{z}}
    ) %>%
    column_to_rownames(rlang::as_name(rlang::ensym(y))) %>%
    as.matrix()
  # browser()
  clustering <- cluster_mat(t(mat))
  df_clustered <- df %>%
    mutate(
      across(
        c({{x}}),
        ~factor(.x, levels = clustering$labels[clustering$order])
      )
    )
  if (is.null(attr(df_clustered, "clusters")))
    attr(df_clustered, "clusters") <- list()
  attr(df_clustered, "clusters")[[rlang::as_name(rlang::ensym(x))]] <- clustering
  df_clustered
}

fgsea_res <- qread(file.path(PATH_PREFIX, "fgsea_res.qs"))
msigdbr_of_interest <- fread(file.path(PATH_PREFIX, "msigdbr_of_interest.csv.gz")) %>%
  as_tibble()
topgo_res <- qread(file.path(PATH_PREFIX, "topgo_res.qs")) %>%
  mutate(
    term_unique = paste(
      str_sub(Term, 1L, 20L),
      GO.ID,
      sep = " "
    )
  )
topgo_bp_genes <- qread(file.path(PATH_PREFIX, "topgo_bp_genes.qs"))

top_go_top_res_directional <- topgo_res %>%
  transmute(
    comparison_type, comparison, comparison_unique, direction,
    term = term_unique, GO.ID,
    p = as.numeric(str_replace(fisher, "< ", "")),
    significant = p < 0.05
  ) %>%
  pivot_wider(
    names_from = direction,
    values_from = c(p, significant)
  ) %>%
  mutate(
    signed_p = case_when(
      significant_up & significant_down ~ case_when(
        p_up < .1 * p_down ~ -log10(p_up),
        p_down < .1 * p_up ~ log10(p_down),
        TRUE ~ NA_real_
      ),
      significant_up ~ -log10(p_up),
      significant_down ~ log10(p_down),
      TRUE ~ .5 * (-log10(p_up) -log10(p_down))
    )
  )

fgsea_gene_stats_w_symbol <- fread(file.path(PATH_PREFIX, "deseq_res_all.csv.gz")) %>%
  as_tibble() %>%
  mutate(
    gene_name = if_else(
      is.na(hgnc_symbol),
      ensembl_gene_id,
      hgnc_symbol
    )
  )

vemu_treated_vs_untreated <- qread(
  file.path(PATH_PREFIX, "deseq_res_vemu_treated_vs_untreated.qs")
) %>%
  mutate(
    plot_output_name = paste0(
      "plot_out_volcano_",
      cell_line,
      "_",
      concentration_fbs
    ),
    across(condition, ~str_replace_all(.x, "[^[:alnum:]]", "_"))
  ) %>%
  arrange(
    cell_line, concentration_fbs
  )

vemu_linear_res <- qread(
  file.path(PATH_PREFIX, "deseq_res_linear_ordinal_by_cell_line.qs")
)

vemu_treated_vs_untreated_long <- vemu_treated_vs_untreated %>%
  mutate(
    condition_print = paste0(
      cell_line,
      " ",
      concentration_fbs
    )
  ) %>%
  unnest(de_res)

gene_meta <- vemu_treated_vs_untreated$de_res[[1]] %>%
  distinct(ensembl_gene_id, hgnc_symbol) %>%
  mutate(
    gene_name = coalesce(hgnc_symbol, ensembl_gene_id)
  ) %>%
  select(ensembl_gene_id, gene_name)

ensembl_gene_id_2_gene_name <- with(
  gene_meta,
  set_names(gene_name, ensembl_gene_id)
)

vemu_treated_vs_untreated_wide <- vemu_treated_vs_untreated %>%
  select(condition, de_res) %>%
  unnest(de_res) %>%
  mutate(
    neglog10padj = -log10(if_else(padj == 0, .5 * min(padj[padj > 0]), padj))
  ) %>%
  select(condition, ensembl_gene_id, hgnc_symbol, log2FoldChange, neglog10padj, padj) %>%
  inner_join(
    gene_meta,
    by = "ensembl_gene_id"
  ) %>%
  mutate(
      hovertext = paste0(
      "<i>Gene:</i> ", gene_name, "<br>",
      "<i>Log2 fold change:</i> ", signif(log2FoldChange, 2), "<br>",
      "<i>Adjusted p-value:</i> ", signif(padj, 2)
    )
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(log2FoldChange, neglog10padj, padj, hovertext)
  )

vemu_linear_res_long <- vemu_linear_res %>%
  select(-n_de) %>%
  unnest(res) %>%
  inner_join(
    gene_meta,
    by = "ensembl_gene_id"
  )

normalized_counts <- fread(
  here("vemu_shiny", "data", "normalized_counts.csv.gz")
) %>%
  as_tibble()

vst_counts <- fread(
  here("vemu_shiny", "data", "counts_vst.csv.gz")
) %>%
  as_tibble()

de_meta_vemu <- fread(
  here("vemu_shiny", "data", "de_meta_vemu.csv")
) %>%
  as_tibble()

filter_by_gene <- function(data, gene_symbol) {
  filtered <- data %>%
    filter(hgnc_symbol == gene_symbol)
  if (nrow(filtered) == 0)
    # Try gene id instead
    filtered <- data %>%
      filter(ensembl_gene_id == gene_symbol)
  if (nrow(filtered) == 0)
    stop("Gene not found")
  filtered
}

plot_de_single_gene_vemu <- function(gene_symbol) {
  filtered <- filter_by_gene(vemu_treated_vs_untreated_long, gene_symbol) %>%
    mutate(
      padj = replace_na(padj, 1),
    )
  p <- ggplot(
    filtered,
    aes(
      concentration_fbs, log2FoldChange,
      color = cell_line,
      shape = padj < 0.05,
      group = cell_line
    )
  ) +
    geom_point() +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 4)
    ) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(shape = "padj < 0.05")
  p
}

plot_expression_single_gene_vemu <- function(gene_symbol, counts = c("normalized", "vst")) {
  counts <- match.arg(counts)
  counts_table <- switch(
    counts,
    vst = vst_counts,
    normalized = normalized_counts
  )
  df <- filter_by_gene(counts_table, gene_symbol)
  p <- df %>%
    pivot_longer(
      -c(ensembl_gene_id, hgnc_symbol), names_to = "well", values_to = "normalized_count"
    ) %>%
    inner_join(
      de_meta_vemu,
      by = "well"
    ) %>%
    mutate(across(starts_with("concentration"), ~fct_inseq(as.character(.x), ordered = TRUE))) %>%
    ggplot(
      aes(concentration, y = normalized_count, fill = concentration_fbs, group = concentration_fbs)
    ) +
      stat_summary(geom = "col", position = "dodge", fun.data = mean_se) +
      geom_quasirandom(dodge.width = 0.9, varwidth = TRUE, position = "dodge", groupOnX = TRUE) +
      facet_wrap(~cell_line) +
      labs(
        x = "Vemurafenib treatment",
        fill = "FBS concentration",
        y = paste(switch(counts, vst = "Variance stabilized", normalized = "Normalized"), "gene count")
      )
   if (counts == "normalized")
     p <- p + scale_y_continuous(
       trans = scales::pseudo_log_trans(base = 10), breaks = function(...) {
         l <- list(...)
         l[[1]][1] <- max(l[[1]][1], 1)
         x <- rlang::exec(scales::breaks_log(), !!!l)
         c(0, x)
       }
      )
   p
}

make_volcano_plot <- function(data, condition, title) {
  plot_ly(
    data,
    x = reformulate(paste0("log2FoldChange_", condition)),
    y = reformulate(paste0("neglog10padj_", condition)),
    customdata = ~ensembl_gene_id,
    hovertext = reformulate(paste0("hovertext_", condition)),
    hoverinfo = "text",
    # hovertemplate = paste(
    #   "<i>Gene:</i> %{gene_name}<br>",
    #   "<i>log2 fold change:</i> %{x}<br>",
    #   "<i>adjusted p-value:</i> %{text}"
    # ),
    type = "scatter",
    mode = "markers",
    marker = list(
      opacity = 0.8
    )
  ) %>%
    layout(
      annotations = list(
        x = 0.5 , y = 1, text = title, showarrow = FALSE,
        xref = 'paper', yref = 'paper',
        xanchor = "center", yanchor = "bottom"
      ),
      showlegend = FALSE
    )
}

ui <- fluidPage(
  actionButton("plot_gene", "Plot selected gene"),
  tabsetPanel(
    type = "tabs",
    tabPanel(
      "Vemurafenib vs. untreated",
      plotlyOutput("plot_out_volcano", height = "600px")
    ),
    tabPanel(
      "Linear model",
      plotlyOutput("plot_out_linear", height = "400px")
    ),
    tabPanel(
      "Selected genes",
      uiOutput("out_selected_genes")
    ),
    tabPanel(
      "Gene set enrichment",
      fluidRow(
        sliderInput(
          "n_significant_per_term",
          "Number of pathways to show per comparison",
          min = 1, max = 100, value = 10
        )
      ),
      fluidRow(
        column(
          width = 6,
          tabsetPanel(
            type = "tabs",
            id = "tabs_gene_set_enrichment",
            tabPanel(
              "FGSEA",
              plotlyOutput("plot_out_fgsea", height = "800px"),
              value = "fgsea"
            ),
            tabPanel(
              "TOPGO",
              plotlyOutput("plot_out_topgo", height = "800px"),
              value = "topgo"
            )
          )
        ),
        column(
          width = 6,
          tabsetPanel(
            type = "tabs",
            tabPanel(
              "Gene set heatmap",
              plotlyOutput("plot_out_gene_set_heatmap", height = "800px")
              # InteractiveComplexHeatmapOutput(
              #   "out_heatmap_gene_expression"
              # )
            ),
            tabPanel(
              "Gene set beeswarm",
              plotOutput("plot_out_gene_set_beeswarm", height = "800px")
            )
          )
        )
      )
    )
  ),
  fluidRow(
    column(
      width = 4,
      uiOutput("out_selected_gene"),
    ),
    column(
      width = 8,
      radioButtons(
        "count_type",
        "Count normalization",
        choices = c(
          "normalized" = "normalized",
          "variance stabilized" = "vst"
        ),
        selected = "normalized"
      )
    )
  ),
  fluidRow(
    column(
      width = 8,
      plotOutput("plot_out_expression_bars")
    ),
    column(
      width = 4,
      plotOutput("plot_out_de")
    )
  )
)

library(crosstalk)
vemu_treated_vs_untreated_ct <- SharedData$new(
  vemu_treated_vs_untreated_wide,
  key = ~gene_name,
  group = "Gene symbol"
)
volcano_plots <- vemu_treated_vs_untreated %>%
  pmap(
    function(cell_line, concentration_fbs, condition, ...) {
      make_volcano_plot(
        data = vemu_treated_vs_untreated_ct,
        condition = condition,
        title = paste0(cell_line, " ", concentration_fbs, "% FBS")
      )
    }
  )

volcano_plots_combined <- subplot(
  volcano_plots,
  nrows = 3,
  shareX = TRUE,
  shareY = FALSE,
  titleX = FALSE,
  titleY = FALSE
) %>%
  layout(
    scene = list(
      xaxis = list(title = "log2FoldChange"),
      yaxis = list(title = "-log10(padj)")
    )
  ) %>%
  highlight(
    on = "plotly_hover",
    persistent = FALSE,
    selectize = list(
      onItemAdd = htmlwidgets::JS(r"{
        function(value, $item) {
          Shiny.setInputValue("hovered_gene", value);
        }
      }"),
      closeAfterSelect = TRUE
    ),
    opacityDim = 1,
    selected = attrs_selected(
      opacity = 1,
      marker = list(color = "black", size = 10)
    ),
    debounce = 10
  ) %>%
  event_register("plotly_selected") %>%
  toWebGL()

volcano_plots_combined$x$source <- "plot_out_volcano"

vemu_linear_res_ct <- vemu_linear_res_long %>%
  select(ensembl_gene_id, gene_name, coef, cell_line, log2FoldChange, padj) %>%
  pivot_wider(names_from = cell_line, values_from = c(log2FoldChange, padj)) %>%
  mutate(
    signif = replace_na(padj_A375 < 0.05, FALSE) | replace_na(padj_WM989 < 0.05, FALSE)
  ) %>%
  arrange(signif) %>%
  SharedData$new(key = ~gene_name, group = "Gene symbol ")

linear_res_plots_combined <- vemu_linear_res_ct %>%
  ggplot(
    aes(
      log2FoldChange_A375, log2FoldChange_WM989,
      alpha = signif,
      color = signif,
      customdata = ensembl_gene_id,
      # for plotly
      text = paste0(
        "<b>", gene_name, "</b><br>",
        "p A375 ", signif(padj_A375, 3), " p WM989 ", signif(padj_WM989, 3), "<br>",
        "log2FC A375 ", signif(log2FoldChange_A375, 3), " log2FC WM989 ", signif(log2FoldChange_WM989, 3)
      )
    )
  ) +
    geom_point(shape = 16) +
    scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.2)) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "gray50")) +
    facet_wrap(~coef, scales = "free")

linear_res_plots_combined_plotly <- linear_res_plots_combined %>%
  ggplotly(source = "plot_out_volcano", tooltip = "text") %>%
  highlight(
    on = "plotly_hover",
    persistent = FALSE,
    selectize = list(
      onItemAdd = htmlwidgets::JS(r"{
        function(value, $item) {
          Shiny.setInputValue("hovered_gene", value);
        }
      }"),
      closeAfterSelect = TRUE
    ),
    opacityDim = 1,
    selected = attrs_selected(
      opacity = 1,
      marker = list(color = "red", size = 10)
    ),
    debounce = 10
  ) %>%
  event_register("plotly_selected") %>%
  toWebGL()

server <- function(input, output, session) {
  r_fgsea_top_pathways <- reactive({
    fgsea_res %>%
      filter(
        pathway %in% {
          msigdbr_of_interest %>%
            # filter(gs_cat == "H") %>%
            pull(gs_name)
        }
      ) %>%
      group_by(comparison, comparison_unique) %>%
      filter(padj < 0.05) %>%
      arrange(pval) %>%
      slice_head(n = input$n_significant_per_term) %>%
      ungroup() %>%
      pull(pathway) %>%
      unique()
  })
  r_fgsea_clustered <- reactive({
    fgsea_res %>%
      filter(pathway %in% r_fgsea_top_pathways()) %>%
      cluster_df(comparison_unique, pathway, NES) %>%
      cluster_df(pathway, comparison_unique, NES)
  })
  r_fgsea_plot <- reactive({
    p <- ggplot(
      r_fgsea_clustered(),
      aes(
        comparison_unique,
        pathway,
        fill = NES,
        customdata = pathway
      )
    ) +
      geom_raster() +
      geom_text(
        aes(label = padj_text),
        data = ~.x %>%
          mutate(
            padj_text = cut(
              padj,
              breaks = c(0, 0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "")
            )
          ),
        size = 2
      ) +
      scale_fill_distiller(palette = "RdBu") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # aspect.ratio = 10
      ) +
      # coord_equal() +
      labs(x = NULL, y = NULL) +
      ggforce::facet_row(~comparison_type, scales = "free_x", space = "free", drop = TRUE)
    ggplotly(p, source = "plotly_fgsea") %>%
      event_register("plotly_click")
  })
  r_topgo_top_pathways <- reactive({
    top_go_top_res_directional %>%
      group_by(comparison, comparison_unique) %>%
      filter(abs(signed_p) > -log10(0.05)) %>%
      arrange(desc(abs(signed_p))) %>%
      slice_head(n = input$n_significant_per_term) %>%
      ungroup() %>%
      pull(term) %>%
      unique()
  })
  r_topgo_clustered <- reactive({
    top_go_top_res_directional %>%
      filter(term %in% r_topgo_top_pathways()) %>%
      cluster_df(comparison_unique, term, signed_p) %>%
      cluster_df(term, comparison_unique, signed_p)
  })
  r_topgo_plot <- reactive({
    plot_range <- find_scale_range(r_topgo_clustered()$signed_p, 0.05)
    # browser()
    p <- ggplot(
      r_topgo_clustered() %>%
        mutate(across(signed_p, ~replace_na(.x, 5))),
      aes(
        comparison_unique,
        term,
        fill = signed_p,
        customdata = GO.ID,
        text = paste0(
          "<b>", term, "<br>", comparison_unique, "</b><br>",
          "p up ", signif(p_up, 3), " p down ", signif(p_down, 3), "<br>",
          "p signed ", signif(signed_p, 3)
        )
      )
    ) +
      geom_raster() +
      scale_fill_distiller(
        palette = "RdBu", type = "div", limits = plot_range, na.value = "gray50",
        oob = scales::squish
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # aspect.ratio = 10
      ) +
      # coord_equal() +
      labs(x = NULL, y = NULL) +
      ggforce::facet_row(~comparison_type, scales = "free_x", space = "free", drop = TRUE)
    ggplotly(p, source = "plotly_topgo", tooltip = "text") %>%
      event_register("plotly_click")
  })
  r_selected_genes_clustered <- reactive({
    req(r_selected_genes_enrichment())
    print(r_selected_genes_enrichment())
    # browser()
    fgsea_gene_stats_w_symbol %>%
      filter(
        ensembl_gene_id %in% r_selected_genes_enrichment()
      ) %>%
      # cluster_df(comparison_unique, gene_symbol, log2FoldChange) %>%
      mutate(
        across(
          comparison_unique,
          \(x) factor(
            x,
            levels = levels(
              isolate(
                if (input$tabs_gene_set_enrichment == "fgsea")
                  r_fgsea_clustered()$comparison_unique
                else
                  r_topgo_clustered()$comparison_unique
              )
            )
          )
        )
      ) %>%
      cluster_df(gene_name, comparison_unique, log2FoldChange)
  })
  r_selected_genes_heatmap <- reactive({
    plot_range <- find_scale_range(r_selected_genes_clustered()$log2FoldChange, 0.05)
    ggplot(
      r_selected_genes_clustered(),
      aes(
        comparison_unique,
        gene_name,
        fill = log2FoldChange
      )
    ) +
      geom_tile() +
      geom_text(
        aes(label = padj_text),
        data = ~.x %>%
          mutate(
            padj_text = cut(
              padj,
              breaks = c(0, 0.001, 0.01, 0.05, 1),
              labels = c("***", "**", "*", "")
            )
          ),
        size = 2
      ) +
      scale_fill_distiller(
        palette = "RdBu", type = "div", limits = plot_range, na.value = "gray50",
        oob = scales::squish
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        # aspect.ratio = 10
      ) +
      # coord_equal() +
      labs(x = NULL, y = NULL) +
      ggforce::facet_row(~comparison_type, scales = "free_x", space = "free", drop = TRUE)
  })
  library(ggbeeswarm)
  r_selected_genes_beeswarm <- reactive({
    ggplot(
      r_selected_genes_clustered() %>%
        mutate(across(comparison_unique, fct_rev)),
      aes(
        log2FoldChange,
        comparison_unique
      )
    ) + 
      geom_quasirandom(
        orientation = "y"
      ) +
      ggforce::facet_col(~comparison_type, scales = "free_y", space = "free", drop = TRUE)
  })
  r_selected_genes_enrichment <- reactiveVal()
  observe({
    new_fgsea_pathway <- event_data("plotly_click", source = "plotly_fgsea")$customdata
    # browser()
    req(new_fgsea_pathway)
    new_gene_set <- msigdbr_of_interest %>%
      filter(gs_name == new_fgsea_pathway) %>%
      pull(ensembl_gene)
    if (length(new_gene_set) > 0)
      r_selected_genes_enrichment(new_gene_set)
  })
  observe({
    click_data <- event_data("plotly_click", source = "plotly_topgo")
    req(click_data)
    # browser()
    # For some reason customdata doesn't come through! Have to use x, y coordinates
    # directly
    new_topgo_pathway <- isolate({
      filter(
        r_topgo_clustered(),
        term == levels(r_topgo_clustered()$term)[click_data$y]
      ) %>%
        chuck("GO.ID", 1)
    })
    print(new_topgo_pathway)
    req(new_topgo_pathway)
    new_gene_set <- topgo_bp_genes[[new_topgo_pathway]]
    if (length(new_gene_set) > 0)
      r_selected_genes_enrichment(new_gene_set)
  })
  output$plot_out_volcano <- renderPlotly(
    volcano_plots_combined
  )
  output$plot_out_linear <- renderPlotly(
    linear_res_plots_combined_plotly
  )
  output$plot_out_fgsea <- renderPlotly(
    r_fgsea_plot()
  )
  output$plot_out_topgo <- renderPlotly(
    r_topgo_plot()
  )
  output$plot_out_gene_set_heatmap <- renderPlotly(
    r_selected_genes_heatmap()
  )
  output$plot_out_gene_set_beeswarm <- renderPlot(
    r_selected_genes_beeswarm()
  )
  r_selected_gene <- reactiveVal()
  observe({
    r_selected_gene(
      ensembl_gene_id_2_gene_name[
        event_data("plotly_click", source = "plot_out_volcano")$customdata
      ]
    )
  })
  bindEvent(
    observe({
      new_gene <- input$hovered_gene
      # Leave reactive unchanged if current selection is invalid
      if(!is.null(new_gene))
        r_selected_gene(new_gene)
    }), input$plot_gene
  )
  r_brushed_genes <- reactive({
    data <- event_data("plotly_selected", source = "plot_out_volcano")
    req(nrow(data) > 0)
    data$customdata
  })
  output$out_selected_genes <- renderUI({
    validate(
      need(r_brushed_genes(), "No genes brushed")
    )
    tagList(
      tags$b("Brushed genes:"),
      rlang::exec(
        tags$ul,
        !!!map(
          r_brushed_genes(),
          ~tags$li(ensembl_gene_id_2_gene_name[.x])
        )
      )
    )
  })
  output$out_selected_gene <- renderUI({
    validate(
      need(r_selected_gene(), "No gene selected")
    )
    tags$b(r_selected_gene())
  })
  output$plot_out_expression_bars <- renderPlot({
    validate(
      need(r_selected_gene(), "Click on a point in the volcano plot to select a gene")
    )
    plot_expression_single_gene_vemu(
      gene_symbol = r_selected_gene(),
      counts = input$count_type
    )
  }, res = 120)
  output$plot_out_de <- renderPlot({
    validate(
      need(r_selected_gene(), "Click on a point in the volcano plot to select a gene")
    )
    plot_de_single_gene_vemu(
      gene_symbol = r_selected_gene()
    )
  }, res = 120)
}

shinyApp(
  ui, server,
  options = list(port = 4321, launch.browser = FALSE)
)

