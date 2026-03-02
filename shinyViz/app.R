library(shiny)
library(Seurat)
library(ggplot2)

# ── Load data once at startup ─────────────────────────────────────────────────
seurat_obj <- readRDS("../dependencies/shared/Post_R3_FINAL_with_counts.rds")

# UMAP coordinates
available_reductions <- names(seurat_obj@reductions)
umap_key <- grep("umap", available_reductions, ignore.case = TRUE, value = TRUE)[1]
if (is.na(umap_key)) stop("No UMAP reduction found in this Seurat object.")

umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = umap_key))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Genes available in scale.data (check all assays)
get_scale_genes <- function(obj) {
  genes <- c()
  for (assay in names(obj@assays)) {
    sd <- GetAssayData(obj, assay = assay, layer = "scale.data")
    if (nrow(sd) > 0) genes <- union(genes, rownames(sd))
  }
  sort(genes)
}
genes <- get_scale_genes(seurat_obj)
if (length(genes) == 0) stop("scale.data layer is empty across all assays.")

# Helper: fetch scaled expression for a gene across assays
get_scaled_expr <- function(obj, gene) {
  for (assay in names(obj@assays)) {
    sd <- GetAssayData(obj, assay = assay, layer = "scale.data")
    if (gene %in% rownames(sd)) return(sd[gene, ])
  }
  return(NULL)
}

# Metadata columns
meta_df   <- seurat_obj@meta.data
meta_cols <- colnames(meta_df)

n_cells <- ncol(seurat_obj)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(

  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Arial, sans-serif; }
    .well { background: #f8f9fa; border: 1px solid #dee2e6; }
    h4 { color: #495057; }
  "))),

  titlePanel("RV Atlas — Single-cell UMAP Viewer"),

  sidebarLayout(
    sidebarPanel(width = 3,

      h4("Color by"),
      radioButtons("color_mode", label = NULL,
                   choices  = c("Gene expression" = "gene",
                                "Metadata"        = "meta"),
                   selected = "gene"),

      conditionalPanel(
        condition = "input.color_mode == 'gene'",
        selectizeInput("gene",
                       label   = "Gene (scale.data):",
                       choices = NULL,
                       options = list(placeholder    = "Type to search…",
                                      maxOptions     = 50,
                                      highlight      = TRUE))
      ),

      conditionalPanel(
        condition = "input.color_mode == 'meta'",
        selectInput("meta_col",
                    label   = "Metadata column:",
                    choices = meta_cols)
      ),

      hr(),
      h4("Display"),
      sliderInput("pt_size", "Point size:",
                  min = 0.05, max = 3, value = 0.4, step = 0.05),
      sliderInput("alpha", "Opacity:",
                  min = 0.1,  max = 1, value = 0.7, step = 0.05),

      conditionalPanel(
        condition = "input.color_mode == 'gene'",
        selectInput("palette", "Color scale:",
                    choices  = c("Blue–White–Red" = "bwr",
                                 "Purple–White–Orange" = "pwo",
                                 "Viridis"        = "viridis"),
                    selected = "bwr")
      ),

      hr(),
      tags$small(paste0("Cells: ", format(n_cells, big.mark = ","),
                        "   |   Genes in scale.data: ",
                        format(length(genes), big.mark = ",")))
    ),

    mainPanel(width = 9,
      plotOutput("umap_plot", height = "680px"),
      downloadButton("dl_plot", "Download PDF")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Server-side selectize so thousands of genes don't bloat the page
  updateSelectizeInput(session, "gene", choices = genes, server = TRUE)

  # Reactive plot object (shared between display and download)
  current_plot <- reactive({

    plot_df <- umap_coords

    if (input$color_mode == "gene") {
      req(input$gene, input$gene %in% genes)

      expr <- get_scaled_expr(seurat_obj, input$gene)
      req(!is.null(expr))
      plot_df$value <- expr[rownames(plot_df)]

      # Sort so highest-expressing cells render on top
      plot_df <- plot_df[order(plot_df$value), ]

      pal <- switch(input$palette,
        bwr     = scale_color_gradient2(low  = "#2166AC",
                                        mid  = "#F7F7F7",
                                        high = "#B2182B",
                                        midpoint = 0,
                                        name = paste0(input$gene, "\n(scaled)")),
        pwo     = scale_color_gradient2(low  = "#762A83",
                                        mid  = "#F7F7F7",
                                        high = "#E08214",
                                        midpoint = 0,
                                        name = paste0(input$gene, "\n(scaled)")),
        viridis = scale_color_viridis_c(name = paste0(input$gene, "\n(scaled)"),
                                        option = "C")
      )

      p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = value)) +
        geom_point(size = input$pt_size, alpha = input$alpha, stroke = 0) +
        pal +
        labs(title = paste("UMAP —", input$gene)) +
        theme_classic(base_size = 13) +
        theme(plot.title    = element_text(face = "bold"),
              axis.line     = element_line(linewidth = 0.4),
              legend.key.height = unit(1.2, "cm"))

    } else {
      req(input$meta_col)
      plot_df$value <- meta_df[rownames(plot_df), input$meta_col]

      p <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = value)) +
        geom_point(size = input$pt_size, alpha = input$alpha, stroke = 0) +
        labs(title = paste("UMAP —", input$meta_col),
             color = input$meta_col) +
        theme_classic(base_size = 13) +
        theme(plot.title = element_text(face = "bold"),
              axis.line  = element_line(linewidth = 0.4))

      # Use discrete palette for low-cardinality columns
      if (is.factor(plot_df$value) || is.character(plot_df$value) ||
          length(unique(plot_df$value)) <= 30) {
        p <- p + scale_color_discrete()
      }
    }

    p
  })

  output$umap_plot <- renderPlot({ current_plot() })

  output$dl_plot <- downloadHandler(
    filename = function() {
      label <- if (input$color_mode == "gene") input$gene else input$meta_col
      paste0("UMAP_", label, ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(),
             width = 9, height = 7, device = "pdf")
    }
  )
}

shinyApp(ui, server)
