library(shiny)
library(plotly)
library(tidyverse)
library(here)
library(qs)
library(ggbeeswarm)

theme_set(theme_minimal())

vemu_treated_vs_untreated <- qread(
  here("data", "deseq_res_vemu_treated_vs_untreated.qs")
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
  here("data", "deseq_res_linear_ordinal_by_cell_line.qs")
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

normalized_counts <- read_csv(
  here("data", "normalized_counts.csv.gz")
)

vst_counts <- read_csv(
  here("data", "counts_vst.csv.gz")
)

de_meta_vemu <- read_csv(
  here("data", "de_meta_vemu.csv")
)

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
  output$plot_out_volcano <- renderPlotly(
    volcano_plots_combined
  )
  output$plot_out_linear <- renderPlotly(
    linear_res_plots_combined_plotly
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
    }),
    input$plot_gene
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

