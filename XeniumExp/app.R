library(shiny)
library(plotly)
library(Seurat)
library(ggplot2)
library(dplyr)

# ── Startup: load & prepare data ─────────────────────────────────────────────
message("Loading Seurat object...")
obj <- readRDS("../dependencies/Figure_3/collated_xenium_data_FINAL.rds")

message("Loading external metadata...")
ext_meta <- read.csv("../dependencies/shared/Xenium_metadata.csv",
                     row.names = 1, check.names = FALSE)

new_cols <- setdiff(colnames(ext_meta), colnames(obj@meta.data))
shared_cells <- intersect(rownames(obj@meta.data), rownames(ext_meta))
obj@meta.data[shared_cells, new_cols] <- ext_meta[shared_cells, new_cols, drop = FALSE]

message("Extracting spatial coordinates...")
fov_images <- grep("^fov", names(obj@images), value = TRUE)

coords_list <- lapply(fov_images, function(img) {
  cc <- tryCatch(
    GetTissueCoordinates(obj, image = img, which = "centroids"),
    error = function(e) NULL
  )
  if (is.null(cc) || nrow(cc) == 0) return(NULL)
  data.frame(
    x       = cc$x,
    y       = cc$y,
    cell_id = cc$cell,
    fov     = img,
    stringsAsFactors = FALSE
  )
})
coords_df <- do.call(rbind, Filter(Negate(is.null), coords_list))
rownames(coords_df) <- coords_df$cell_id

keep_meta <- c("patient", "group", "names", "subnames", "markernames",
               "niche_manual", "niches.01", "niches.05", "niches",
               "niche_snn_res.0.1", "broad_niches.01",
               "nCount_Xenium", "nFeature_Xenium",
               "SCT_snn_res.0.3", "SCT_snn_res.0.5", "SCT_snn_res.1",
               new_cols[grepl("^kmeans_", new_cols)])
keep_meta <- intersect(keep_meta, colnames(obj@meta.data))

present_cells <- intersect(coords_df$cell_id, rownames(obj@meta.data))
coords_df[present_cells, keep_meta] <- obj@meta.data[present_cells, keep_meta]

# Patient selector
patient_info <- unique(coords_df[, c("fov", "patient", "group")])
patient_info <- patient_info[order(patient_info$patient), ]
patient_choices <- setNames(
  patient_info$fov,
  paste0("Patient ", patient_info$patient, "  [", patient_info$group, "]")
)

message("Extracting scale.data...")
scale_mat <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
gene_list  <- sort(rownames(scale_mat))

meta_display_cols <- c(
  "names", "subnames", "markernames", "group", "patient",
  "niche_manual", "niches.01", "niches.05", "niches",
  "niche_snn_res.0.1", "broad_niches.01",
  "nCount_Xenium", "nFeature_Xenium",
  "SCT_snn_res.0.3", "SCT_snn_res.0.5", "SCT_snn_res.1",
  new_cols[grepl("^kmeans_", new_cols)]
)
meta_display_cols <- intersect(meta_display_cols, colnames(coords_df))

# ── Colour helpers ────────────────────────────────────────────────────────────
# 30-colour discrete palette (scCustomize-inspired)
DISCRETE_PALETTE <- c(
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
  "#A65628","#F781BF","#999999","#66C2A5","#FC8D62",
  "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
  "#B3B3B3","#1B9E77","#D95F02","#7570B3","#E7298A",
  "#66A61E","#E6AB02","#A6761D","#666666","#2166AC",
  "#D1E5F0","#92C5DE","#4393C3","#F4A582","#D6604D"
)

discrete_colors <- function(vals) {
  lvls <- unique(na.omit(as.character(vals)))
  lvls <- sort(lvls)
  pal  <- rep(DISCRETE_PALETTE, length.out = length(lvls))
  setNames(pal, lvls)
}

bwr_colorscale <- list(
  list(0,   "#2166AC"),
  list(0.5, "#F7F7F7"),
  list(1,   "#B2182B")
)
pwo_colorscale <- list(
  list(0,   "#762A83"),
  list(0.5, "#F7F7F7"),
  list(1,   "#E08214")
)
viridis_colorscale <- list(
  list(0,    "#440154"),
  list(0.25, "#31688E"),
  list(0.5,  "#35B779"),
  list(0.75, "#FDE725"),
  list(1,    "#FDE725")
)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- tagList(
  tags$head(
    tags$link(rel = "stylesheet", href = "style.css"),
    tags$title("XeniumExp — RV Atlas Spatial Viewer")
  ),

  # ── Header ─────────────────────────────────────────────────────────────────
  tags$div(id = "xe-header",

    # Brand
    tags$div(class = "xe-brand",
      tags$span(class = "xe-accent", "Xenium"),
      tags$span("Exp")
    ),

    # Patient selector (centre)
    tags$div(class = "xe-patient-selector",
      tags$label(`for` = "patient_fov", "Patient"),
      selectInput("patient_fov", label = NULL,
                  choices  = patient_choices,
                  selected = patient_choices[1],
                  width    = "260px")
    ),

    # Header action buttons
    tags$div(class = "xe-header-actions",
      actionButton("info_btn",  "ℹ Info",   class = "btn-header"),
      actionButton("tool_lasso","⌀ Lasso",  class = "btn-header active", id = "tool_lasso"),
      actionButton("tool_box",  "▭ Box",    class = "btn-header",        id = "tool_box")
    )
  ),

  # ── Body ───────────────────────────────────────────────────────────────────
  tags$div(id = "xe-body",

    tags$div(id = "xe-workspace",

      # ── Sidebar ──────────────────────────────────────────────────────────
      tags$div(id = "xe-sidebar",

        tags$div(class = "xe-section-header", "Color By"),

        radioButtons("color_mode", label = NULL,
                     choices  = c("Metadata" = "meta", "Gene expression" = "gene"),
                     selected = "meta"),

        conditionalPanel(
          condition = "input.color_mode == 'meta'",
          selectInput("meta_col", "Metadata column:",
                      choices  = meta_display_cols,
                      selected = "names",
                      width    = "100%")
        ),

        conditionalPanel(
          condition = "input.color_mode == 'gene'",
          selectizeInput("gene", "Gene (SCT scale.data):",
                         choices = NULL,
                         options = list(placeholder = "Type to search…",
                                        maxOptions  = 50),
                         width   = "100%"),
          selectInput("palette", "Color scale:",
                      choices  = c("Blue–White–Red" = "bwr",
                                   "Purple–White–Orange" = "pwo",
                                   "Viridis" = "viridis"),
                      selected = "bwr",
                      width    = "100%")
        ),

        tags$div(class = "xe-section-header", "Display"),

        sliderInput("pt_size", "Point size:",
                    min = 1, max = 8, value = 3, step = 0.5),
        sliderInput("alpha", "Opacity:",
                    min = 0.1, max = 1, value = 0.7, step = 0.05),

        numericInput("max_cells", "Max cells (0 = all):",
                     value = 60000, min = 0, step = 10000,
                     width = "100%"),

        tags$div(class = "xe-section-header", "Selection"),
        actionButton("clear_sel", "Clear Selection",
                     style = "width:100%;font-size:12px;margin-bottom:6px;"),

        uiOutput("cell_count_ui"),

        tags$div(class = "xe-section-header", "Export"),
        downloadButton("dl_plot", "Download PDF",
                       class = "btn-download")
      ),

      # ── Main (plot + stats) ───────────────────────────────────────────────
      tags$div(id = "xe-main",

        # Spatial plot (fills remaining space)
        tags$div(id = "xe-plot-container",
          plotlyOutput("spatial_plot",
                       width  = "100%",
                       height = "100%")
        ),

        # Stats strip
        tags$div(id = "xe-stats",

          tags$div(class = "xe-stats-panel",
            tags$div(class = "xe-stats-title", "Cell Composition"),
            plotlyOutput("composition_plot", width = "100%", height = "210px")
          ),

          tags$div(class = "xe-stats-panel",
            tags$div(class = "xe-stats-title", "Gene Expression"),
            plotlyOutput("violin_plot", width = "100%", height = "210px")
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Server-side selectize for gene list
  updateSelectizeInput(session, "gene", choices = gene_list, server = TRUE)

  # ── Drag mode state ────────────────────────────────────────────────────────
  drag_mode <- reactiveVal("lasso")

  observeEvent(input$tool_lasso, {
    drag_mode("lasso")
    plotlyProxy("spatial_plot", session) %>%
      plotlyProxyInvoke("relayout", list(dragmode = "lasso"))
  })

  observeEvent(input$tool_box, {
    drag_mode("select")
    plotlyProxy("spatial_plot", session) %>%
      plotlyProxyInvoke("relayout", list(dragmode = "select"))
  })

  # ── Clear selection ────────────────────────────────────────────────────────
  observeEvent(input$clear_sel, {
    plotlyProxy("spatial_plot", session) %>%
      plotlyProxyInvoke("relayout", list(selections = list()))
    session$sendCustomMessage("clearPlotlySelection", list())
  })

  # ── Patient cells ──────────────────────────────────────────────────────────
  patient_cells <- reactive({
    req(input$patient_fov)
    coords_df[coords_df$fov == input$patient_fov, , drop = FALSE]
  })

  # ── Plot data (subsampled if needed) ──────────────────────────────────────
  plot_df <- reactive({
    df    <- patient_cells()
    max_c <- input$max_cells
    if (max_c > 0 && nrow(df) > max_c) {
      df <- df[sample(nrow(df), max_c), ]
    }
    df
  })

  # ── Cell count UI ──────────────────────────────────────────────────────────
  output$cell_count_ui <- renderUI({
    n_total    <- nrow(patient_cells())
    n_shown    <- nrow(plot_df())
    sel_ids    <- selected_ids()
    n_sel      <- if (is.null(sel_ids)) 0L else length(sel_ids)
    n_sel_show <- if (is.null(sel_ids)) "" else
      tags$span(class = "xe-count-selected",
                paste0(" | ", format(n_sel, big.mark = ","), " selected"))

    tags$div(class = "xe-cell-count",
      paste0(format(n_total, big.mark = ","), " cells"),
      if (n_shown < n_total)
        paste0(" (", format(n_shown, big.mark = ","), " shown)")
      else "",
      n_sel_show
    )
  })

  # ── Selected cell IDs (from plotly lasso/box) ─────────────────────────────
  selected_ids <- reactive({
    d <- event_data("plotly_selected", source = "spatial")
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d$key
  })

  # ── Spatial plot ──────────────────────────────────────────────────────────
  output$spatial_plot <- renderPlotly({
    df <- plot_df()
    req(nrow(df) > 0)

    pt  <- unique(df$patient)[1]
    grp <- unique(df$group)[1]

    # ── Build colour mapping ──────────────────────────────────────────────
    if (input$color_mode == "gene") {
      req(input$gene, input$gene != "", input$gene %in% gene_list)
      expr <- as.numeric(scale_mat[input$gene, df$cell_id])

      # Rescale to [0, 1] for plotly colorscale
      expr_min <- min(expr, na.rm = TRUE)
      expr_max <- max(expr, na.rm = TRUE)
      expr_norm <- if (expr_max > expr_min) {
        (expr - expr_min) / (expr_max - expr_min)
      } else {
        rep(0.5, length(expr))
      }

      cs <- switch(input$palette,
        bwr     = bwr_colorscale,
        pwo     = pwo_colorscale,
        viridis = viridis_colorscale
      )

      # Order: low expression behind, high on top
      ord <- order(expr)
      df  <- df[ord, ]
      expr_norm <- expr_norm[ord]
      raw_expr  <- expr[ord]

      tooltip <- paste0(
        "<b>", df$cell_id, "</b><br>",
        "Type: ", df$names, "<br>",
        input$gene, ": ", round(raw_expr, 3)
      )

      p <- plot_ly(
        x    = ~df$x,
        y    = ~df$y,
        type = "scattergl",
        mode = "markers",
        marker = list(
          color   = expr_norm,
          colorscale = cs,
          size    = input$pt_size,
          opacity = input$alpha,
          line    = list(width = 0),
          colorbar = list(
            title      = paste0(input$gene, "<br>(scaled)"),
            thickness  = 14,
            len        = 0.6,
            tickfont   = list(size = 10)
          ),
          showscale = TRUE
        ),
        text      = tooltip,
        hoverinfo = "text",
        key       = ~df$cell_id,
        source    = "spatial"
      )

    } else {
      # Metadata / discrete coloring
      req(input$meta_col, input$meta_col %in% colnames(df))
      col_vals  <- as.character(df[[input$meta_col]])
      pal_map   <- discrete_colors(col_vals)
      dot_colors <- pal_map[col_vals]
      dot_colors[is.na(dot_colors)] <- "#cccccc"

      tooltip <- paste0(
        "<b>", df$cell_id, "</b><br>",
        input$meta_col, ": ", col_vals
      )

      # Build one trace per category for a proper legend
      cats <- names(pal_map)
      traces <- lapply(cats, function(cat) {
        idx <- which(col_vals == cat)
        if (length(idx) == 0) return(NULL)
        list(
          x         = df$x[idx],
          y         = df$y[idx],
          name      = cat,
          marker_color = pal_map[[cat]],
          text      = tooltip[idx],
          cell_ids  = df$cell_id[idx]
        )
      })
      traces <- Filter(Negate(is.null), traces)

      p <- plot_ly(source = "spatial")
      for (tr in traces) {
        p <- add_trace(p,
          x    = tr$x,
          y    = tr$y,
          type = "scattergl",
          mode = "markers",
          name = tr$name,
          marker = list(
            color   = tr$marker_color,
            size    = input$pt_size,
            opacity = input$alpha,
            line    = list(width = 0)
          ),
          text      = tr$text,
          hoverinfo = "text",
          key       = tr$cell_ids,
          showlegend = TRUE
        )
      }
    }

    p %>% layout(
      title  = list(
        text   = paste0("Patient ", pt, "  |  Group: ", grp,
                        "  |  ", input$color_mode),
        font   = list(size = 13, color = "#374151"),
        x      = 0.02, xanchor = "left"
      ),
      dragmode   = drag_mode(),
      xaxis = list(
        scaleanchor = "y", scaleratio = 1,
        title = "x (µm)", tickfont = list(size = 10),
        zeroline = FALSE, gridcolor = "#e5e7eb"
      ),
      yaxis = list(
        autorange = "reversed",
        title = "y (µm)", tickfont = list(size = 10),
        zeroline = FALSE, gridcolor = "#e5e7eb"
      ),
      showlegend     = (input$color_mode == "meta"),
      legend = list(
        font         = list(size = 10),
        itemsizing   = "constant",
        tracegroupgap = 2
      ),
      paper_bgcolor = "#f8f9fa",
      plot_bgcolor  = "#f8f9fa",
      margin        = list(l = 55, r = 10, t = 35, b = 45)
    ) %>%
      event_register("plotly_selected") %>%
      config(
        scrollZoom   = TRUE,
        displaylogo  = FALSE,
        modeBarButtonsToRemove = c("sendDataToCloud", "resetScale2d",
                                   "hoverClosestCartesian", "hoverCompareCartesian",
                                   "toggleSpikelines")
      )
  })

  # ── Composition plot ───────────────────────────────────────────────────────
  output$composition_plot <- renderPlotly({
    df      <- plot_df()
    sel_ids <- selected_ids()
    req(nrow(df) > 0, "names" %in% colnames(df))

    cell_types <- as.character(df$names)
    pal_map    <- discrete_colors(cell_types)
    all_lvls   <- names(pal_map)

    all_counts  <- table(factor(cell_types, levels = all_lvls))
    all_pct     <- round(100 * all_counts / sum(all_counts), 1)

    if (!is.null(sel_ids)) {
      sel_in_view  <- intersect(sel_ids, df$cell_id)
      sel_types    <- as.character(df$names[df$cell_id %in% sel_in_view])
      sel_counts   <- table(factor(sel_types, levels = all_lvls))
      sel_pct      <- round(100 * sel_counts / sum(sel_counts + 1e-9), 1)

      p <- plot_ly()
      p <- add_trace(p,
        x          = as.numeric(all_pct),
        y          = all_lvls,
        type       = "bar",
        orientation = "h",
        name       = "All cells",
        marker     = list(color = unname(pal_map), opacity = 0.4),
        hovertemplate = "%{y}: %{x}%<extra>All</extra>"
      )
      p <- add_trace(p,
        x          = as.numeric(sel_pct),
        y          = all_lvls,
        type       = "bar",
        orientation = "h",
        name       = "Selected",
        marker     = list(color = unname(pal_map), opacity = 1),
        hovertemplate = "%{y}: %{x}%<extra>Selected</extra>"
      )
      p <- layout(p, barmode = "overlay")
    } else {
      p <- plot_ly(
        x          = as.numeric(all_pct),
        y          = all_lvls,
        type       = "bar",
        orientation = "h",
        name       = "All cells",
        marker     = list(color = unname(pal_map), opacity = 0.85),
        hovertemplate = "%{y}: %{x}%<extra></extra>"
      )
    }

    p %>% layout(
      xaxis = list(title = "%", tickfont = list(size = 9), zeroline = FALSE),
      yaxis = list(title = "", tickfont = list(size = 9), autorange = "reversed"),
      showlegend   = !is.null(sel_ids),
      legend       = list(font = list(size = 9), orientation = "h", y = -0.15),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin        = list(l = 10, r = 10, t = 5, b = 30)
    ) %>%
      config(displayModeBar = FALSE)
  })

  # ── Violin / expression plot ───────────────────────────────────────────────
  output$violin_plot <- renderPlotly({
    if (input$color_mode != "gene" || is.null(input$gene) || input$gene == "") {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(
            annotations = list(list(
              text      = "Select a gene to see expression distribution",
              showarrow = FALSE,
              font      = list(size = 13, color = "#9ca3af"),
              xref = "paper", yref = "paper", x = 0.5, y = 0.5
            )),
            paper_bgcolor = "white", plot_bgcolor = "white"
          ) %>%
          config(displayModeBar = FALSE)
      )
    }

    req(input$gene %in% gene_list)

    df      <- plot_df()
    sel_ids <- selected_ids()
    req(nrow(df) > 0, "names" %in% colnames(df))

    expr       <- as.numeric(scale_mat[input$gene, df$cell_id])
    cell_types <- as.character(df$names)
    all_lvls   <- sort(unique(cell_types[!is.na(cell_types)]))

    p <- plot_ly()

    # All cells (grey violin)
    for (ct in all_lvls) {
      idx <- which(cell_types == ct)
      if (length(idx) < 3) next
      p <- add_trace(p,
        y          = expr[idx],
        type       = "violin",
        name       = ct,
        legendgroup = ct,
        showlegend  = FALSE,
        box        = list(visible = FALSE),
        meanline   = list(visible = FALSE),
        line       = list(color = "#9ca3af", width = 1),
        fillcolor  = "#d1d5db",
        opacity    = 0.6,
        x0         = ct,
        hovertemplate = paste0(ct, "<br>%{y:.2f}<extra>All</extra>")
      )
    }

    # Selected cells (orange)
    if (!is.null(sel_ids)) {
      sel_in_view <- intersect(sel_ids, df$cell_id)
      sel_types   <- cell_types[df$cell_id %in% sel_in_view]
      sel_expr    <- expr[df$cell_id %in% sel_in_view]

      for (ct in all_lvls) {
        idx <- which(sel_types == ct)
        if (length(idx) < 3) next
        p <- add_trace(p,
          y          = sel_expr[idx],
          type       = "violin",
          name       = paste0(ct, " (sel)"),
          legendgroup = paste0(ct, "_sel"),
          showlegend  = FALSE,
          box        = list(visible = FALSE),
          meanline   = list(visible = TRUE, color = "#ea580c", width = 2),
          line       = list(color = "#f97316", width = 1.5),
          fillcolor  = "#fed7aa",
          opacity    = 0.8,
          x0         = ct,
          hovertemplate = paste0(ct, "<br>%{y:.2f}<extra>Selected</extra>")
        )
      }
    }

    p %>% layout(
      violinmode = "overlay",
      xaxis = list(title = "",    tickfont = list(size = 9),
                   tickangle = -30, zeroline = FALSE),
      yaxis = list(title = paste0(input$gene, " (scaled)"),
                   tickfont = list(size = 9), zeroline = TRUE,
                   zerolinecolor = "#e5e7eb"),
      showlegend    = FALSE,
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin        = list(l = 50, r = 10, t = 5, b = 45)
    ) %>%
      config(displayModeBar = FALSE)
  })

  # ── Info modal ────────────────────────────────────────────────────────────
  observeEvent(input$info_btn, {
    df  <- patient_cells()
    pt  <- unique(df$patient)[1]
    grp <- unique(df$group)[1]
    fov <- input$patient_fov
    n   <- nrow(df)

    showModal(modalDialog(
      title = tags$span(style = "color:white;font-weight:700;", "Sample Information"),
      tags$table(class = "info-table",
        tags$tr(tags$td("Patient"),  tags$td(pt)),
        tags$tr(tags$td("Group"),    tags$td(grp)),
        tags$tr(tags$td("FOV"),      tags$td(fov)),
        tags$tr(tags$td("Cells"),    tags$td(format(n, big.mark = ","))),
        tags$tr(tags$td("Genes"),    tags$td(format(length(gene_list), big.mark = ","))),
        tags$tr(tags$td("Assay"),    tags$td("SCT scale.data"))
      ),
      footer = modalButton("Close"),
      easyClose = TRUE,
      size = "s"
    ))
  })

  # ── Download ──────────────────────────────────────────────────────────────
  output$dl_plot <- downloadHandler(
    filename = function() {
      df    <- patient_cells()
      pt    <- unique(df$patient)[1]
      label <- if (input$color_mode == "gene") input$gene else input$meta_col
      paste0("Xenium_P", pt, "_", label, ".pdf")
    },
    content = function(file) {
      df <- plot_df()
      req(nrow(df) > 0)

      if (input$color_mode == "gene") {
        req(input$gene, input$gene %in% gene_list)
        df$color_val <- as.numeric(scale_mat[input$gene, df$cell_id])
        df <- df[order(df$color_val), ]
        cs <- switch(input$palette,
          bwr     = scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                          high = "#B2182B", midpoint = 0,
                                          name = paste0(input$gene, "\n(scaled)")),
          pwo     = scale_color_gradient2(low = "#762A83", mid = "#F7F7F7",
                                          high = "#E08214", midpoint = 0,
                                          name = paste0(input$gene, "\n(scaled)")),
          viridis = scale_color_viridis_c(name = paste0(input$gene, "\n(scaled)"),
                                          option = "C")
        )
        p <- ggplot(df, aes(x, y, color = color_val)) +
          geom_point(size = 0.2, alpha = input$alpha, stroke = 0) +
          cs + scale_y_reverse() + coord_equal()

      } else {
        req(input$meta_col %in% colnames(df))
        df$color_val <- as.character(df[[input$meta_col]])
        pal_map      <- discrete_colors(df$color_val)
        p <- ggplot(df, aes(x, y, color = color_val)) +
          geom_point(size = 0.2, alpha = input$alpha, stroke = 0) +
          scale_color_manual(values = pal_map, name = input$meta_col, na.value = "grey85") +
          scale_y_reverse() + coord_equal()
      }

      pt  <- unique(df$patient)[1]
      grp <- unique(df$group)[1]
      p <- p + theme_void(base_size = 12) +
        labs(subtitle = paste0("Patient ", pt, "  |  Group: ", grp)) +
        theme(
          plot.subtitle     = element_text(color = "grey40", margin = margin(b = 8)),
          legend.position   = "right",
          legend.title      = element_text(size = 10),
          legend.text       = element_text(size = 8),
          plot.background   = element_rect(fill = "white", color = NA),
          plot.margin       = margin(10, 10, 10, 10)
        )

      ggsave(file, plot = p, width = 11, height = 9, device = "pdf")
    }
  )
}

shinyApp(ui, server)
