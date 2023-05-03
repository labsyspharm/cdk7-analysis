library(shiny)
library(bslib)
library(plotly)
library(tidyverse)
library(here)
library(qs)
library(ggbeeswarm)

vemu_treated_vs_untreated <- qread(
  here("deseq", "deseq_res_vemu_treated_vs_untreated.qs")
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
  select(condition, ensembl_gene_id, hgnc_symbol, log2FoldChange, neglog10padj) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(log2FoldChange, neglog10padj)
  ) %>%
  inner_join(
    gene_meta,
    by = "ensembl_gene_id"
  )

normalized_counts <- read_csv(
  here("qc", "normalized_counts.csv.gz")
)

vst_counts <- read_csv(
  here("qc", "counts_vst.csv.gz")
)

de_meta_vemu <- read_csv(
  here("deseq", "de_meta_vemu.csv")
)

plot_expression_single_gene_vemu <- function(gene_symbol, counts = c("normalized", "vst")) {
  counts <- match.arg(counts)
  counts_table <- switch(
    counts,
    vst = vst_counts,
    normalized = normalized_counts
  )
  df <- counts_table %>%
    filter(hgnc_symbol == gene_symbol)
  if (nrow(df) == 0)
    # Try gene id instead
    df <- counts_table %>%
      filter(ensembl_gene_id == gene_symbol)
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
      geom_quasirandom(dodge.width = 0.9, varwidth = TRUE, position = "dodge") +
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
    hovertext = ~gene_name,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    # make the points a bit transparent
    opacity = 0.5
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
  plotlyOutput("plot_out_volcano", height = "600px"),
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
  plotOutput("plot_out_expression_bars")
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
  shareY = TRUE,
  titleX = FALSE,
  titleY = FALSE
) %>%
  layout(
    scene = list(
      xaxis = list(title = "log2FoldChange"),
      yaxis = list(title = "-log10(padj)")
    )
  ) %>%
  event_register("plotly_click") %>%
  highlight(
    on = "plotly_hover",
    persistent = FALSE,
    selectize = TRUE,
    selected = attrs_selected(
      opacity = 1,
      marker = list(color = "black", size = 10)
    )
  ) %>%
  toWebGL()

volcano_plots_combined$x$source <- "plot_out_volcano"

server <- function(input, output, session) {
  output$plot_out_volcano <- renderPlotly(
    volcano_plots_combined
  )
  r_selected_gene <- reactive({
    event_data("plotly_click", source = "plot_out_volcano")$customdata
  })
  output$out_selected_gene <- renderUI({
    req(r_selected_gene())
    tags$b(ensembl_gene_id_2_gene_name[r_selected_gene()])
  })
  output$plot_out_expression_bars <- renderPlot({
    validate(
      need(r_selected_gene(), "Click on a point in the volcano plot to select a gene")
    )
    plot_expression_single_gene_vemu(
      gene_symbol = r_selected_gene(),
      counts = input$count_type
    )
  })
}

shinyApp(
  ui, server,
  options = list(port = 4321, launch.browser = FALSE)
)






