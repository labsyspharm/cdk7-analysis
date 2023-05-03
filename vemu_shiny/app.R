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
    )
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
    #  {
    #    if (counts == "normalized")
    #      mutate(
    #        .,
    #        zero_count = normalized_count == 0,
    #        normalized_count = if_else(zero_count, counts_normalized_nzmin, normalized_count)
    #       )
    #    else mutate(., zero_count = FALSE)
    # } %>%
    # group_by(cell_line, agent, concentration, concentration_fbs, timepoint) %>%
    # summarize(across(normalized_count, mean), .groups = "drop") %>%
    mutate(across(starts_with("concentration"), ~fct_inseq(as.character(.x), ordered = TRUE))) %>%
    ggplot(
      aes(concentration, y = normalized_count, fill = concentration_fbs, group = concentration_fbs)
    ) +
      stat_summary(geom = "col", position = "dodge") +
      geom_quasirandom(dodge.width = 0.9, varwidth = TRUE, position = "dodge") +
      facet_wrap(~cell_line) +
      # scale_shape_manual(values = c(`TRUE` = 25, `FALSE` = 16), guide = "none") +
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
     # p <- p + scale_y_log10()
   p
}



make_volcano_plot <- function(res) {
  plot_ly(
    res %>%
      drop_na(log2FoldChange, padj),
    x = ~log2FoldChange,
    y = ~-log10(padj),
    customdata = ~hgnc_symbol,
    hovertext = ~hgnc_symbol,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    width = 200, height = 200
  ) %>%
    htmlwidgets::onRender(
      "
      function(el) {
        el.on('plotly_hover', function(d) {
          var id = d.points[0].customdata;
          console.log(id);
          PlotlyHoverGroups.setHoveredPoint('volcano-plots', id);
        });
      }
    ") %>%
    toWebGL()
}

        # el.on('plotly_unhover', function(d) {
        #   PlotlyHoverGroups.resetHoveredPoints('volcano-plots');
        # });

volcano_plot_outputs <- pmap(
  vemu_treated_vs_untreated,
  function(plot_output_name, cell_line, concentration_fbs, de_res, ...) {
    card(
      card_header(cell_line),
      card_body(
        p(concentration_fbs, "% FBS"),
        plotlyOutput(
          plot_output_name
          # width = "50rem", height = "50rem"
        )
      )
    )
  }
)

ui <- page_fluid(
  theme = bs_theme(version = 5),
  tags$head(tags$script(src="js/plotly_hover_groups.js")),
  verticalLayout(
    rlang::exec(layout_column_wrap, width = "200px", !!!volcano_plot_outputs),
    plotOutput("plot_out_expression_bars", width = "50rem", height = "50rem"),
  )
)

server <- function(input, output, session) {
  pwalk(
    vemu_treated_vs_untreated,
    function(plot_output_name, de_res, ...) {
      output[[plot_output_name]] <- renderPlotly(
        make_volcano_plot(de_res)
      )
    }
  )
  session$sendCustomMessage(
    type = "add-plotly-hover-group",
    message = list(
      id = "volcano-plots",
      members = vemu_treated_vs_untreated$plot_output_name
    )
  )
}

addResourcePath(
  prefix = "js",
  directoryPath = here("vemu_shiny", "js")
)
shinyApp(
  ui, server,
  options = list(port = 4321, launch.browser = FALSE)
)



