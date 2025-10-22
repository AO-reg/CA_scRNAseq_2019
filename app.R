install.packages("shiny")
install.packages("SeuratObject")
install.packages("spatstat.data", repos = "https://cloud.r-project.org", type = "source")
install.packages("Seurat")
install.packages("ggplot2")
install.packages("DT")
install.packages("future")


# ==== 必要パッケージ ====
library(shiny)
library(Seurat)
library(DT)
library(future)
library(grid)  # unit() で使用

# ==== 既定の相対パス（同梱 RDS）====
DEFAULT_RDS_PATH <- file.path("EMTAB7716.umap_cosine_walktrap.rds")

# ==== groupby の候補（存在する列のみがUIに出ます）====
ALLOWED_GROUPBY <- c("cluster", "DevelopmentalStage", "DaysPostAmputation", "CellCyclePhase")

# ================= UI =================
ui <- fluidPage(
  titlePanel("scRNA-seq 可視化ツール（Seurat RDS 読み込み）"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 fileInput("rds", "Seuratオブジェクト(.rds)をアップロード", accept = ".rds"),
                 textInput("path", "ローカル/サーバ上のRDSパス（任意）", value = ""),
                 actionButton("loadBtn", "読み込む"),
                 hr(),
                 selectInput("groupby", "グループ化（Identsに設定）", choices = NULL),
                 checkboxInput("setIdents", "Identsを上で選んだ列に設定する", value = TRUE),
                 textInput("genes", "遺伝子名（カンマ区切り）", value = ""),
                 selectInput("reduction", "Reduction（次元縮約）", choices = c("umap","tsne","pca"), selected = "umap"),
                 hr(),
                 numericInput("workers", "並列ワーカー数（コア数）", value = 1, min = 1, step = 1),
                 numericInput("memGB", "メモリ上限（GB）", value = 1, min = 1, step = 1),
                 hr(),
                 checkboxInput("onlypos", "マーカーは正方向のみ（only.pos）", value = TRUE),
                 actionButton("runMarkers", "FindAllMarkers 実行")
    ),
    mainPanel(width = 9,
              tabsetPanel(
                tabPanel("UMAP",        plotOutput("p_dim", height = 800)),
                tabPanel("FeaturePlot", plotOutput("p_feat", height = 550)),
                tabPanel("Violin",      plotOutput("p_vln",  height = 550)),
                tabPanel("DotPlot",     plotOutput("p_dot",  height = 550)),
                tabPanel("Markers",     DTOutput("tbl_markers"))
              )
    )
  )
)

# ================= Server =================
server <- function(input, output, session) {
  
  obj <- reactiveVal(NULL)
  
  # --- 起動直後：同梱RDS（data/object.rds）を自動ロード ---
  session$onFlushed(function(){
    # ★ここを修正：reactive を読むときは isolate() を使う
    if (!is.null(isolate(obj()))) return(NULL)
    if (!file.exists(DEFAULT_RDS_PATH)) return(NULL)
    
    x <- try(readRDS(DEFAULT_RDS_PATH), silent = TRUE)
    if (!inherits(x, "try-error") && inherits(x, "Seurat")) {
      obj(x)
      # （以下はそのまま）
      md_cols <- colnames(x@meta.data)
      choices <- intersect(ALLOWED_GROUPBY, md_cols)
      if ("cluster" %in% md_cols)        choices <- c(choices, "cluster")
      if ("seurat_clusters" %in% md_cols) choices <- c(choices, "seurat_clusters")
      choices <- unique(c(choices, md_cols))
      if (length(choices) == 0) choices <- "orig.ident"
      sel <- if ("cluster" %in% choices) "cluster" else if ("seurat_clusters" %in% choices) "seurat_clusters" else choices[1]
      updateSelectInput(session, "groupby", choices = choices, selected = sel)
      
      red_ok <- intersect(c("umap","tsne","pca"), Reductions(x))
      if (length(red_ok) == 0) red_ok <- "pca"
      updateSelectInput(session, "reduction", choices = red_ok, selected = red_ok[1])
    }
  }, once = TRUE)
  
  # --- 手動読み込み（アップロード or 任意パス） ---
  observeEvent(input$loadBtn, {
    x <- NULL
    if (!is.null(input$rds)) {
      x <- try(readRDS(input$rds$datapath), silent = TRUE)
    }
    if (is.null(x) || inherits(x, "try-error")) {
      p <- trimws(input$path)
      if (nzchar(p)) x <- try(readRDS(p), silent = TRUE)
    }
    validate(need(!is.null(x) && inherits(x, "Seurat"), "Seuratオブジェクトを読み込めませんでした"))
    obj(x)
    
    md_cols <- colnames(x@meta.data)
    choices <- intersect(ALLOWED_GROUPBY, md_cols)
    if ("cluster" %in% md_cols)        choices <- c(choices, "cluster")
    if ("seurat_clusters" %in% md_cols) choices <- c(choices, "seurat_clusters")
    choices <- unique(c(choices, md_cols))
    if (length(choices) == 0) choices <- "orig.ident"
    sel <- if ("cluster" %in% choices) "cluster" else if ("seurat_clusters" %in% choices) "seurat_clusters" else choices[1]
    updateSelectInput(session, "groupby", choices = choices, selected = sel)
    
    red_ok <- intersect(c("umap","tsne","pca"), Reductions(x))
    if (length(red_ok) == 0) red_ok <- "pca"
    updateSelectInput(session, "reduction", choices = red_ok, selected = red_ok[1])
  })
  
  # --- 並列・メモリ設定 ---
  observe({
    plan(multisession, workers = as.integer(input$workers))
    options(future.globals.maxSize = as.numeric(input$memGB) * 1024^3)
  })
  
  # --- Idents の設定（必要に応じて） ---
  current_obj <- reactive({
    x <- obj(); req(x)
    if (isTRUE(input$setIdents) && !is.null(input$groupby) && input$groupby %in% colnames(x@meta.data)) {
      try({ Idents(x) <- x@meta.data[[input$groupby]] }, silent = TRUE)
    }
    x
  })
  
  # --- 遺伝子のパース ---
  genes <- reactive({
    g <- trimws(unlist(strsplit(input$genes, ",")))
    g[g != ""]
  })
  
  # --- UMAP 図 ---
  output$p_dim <- renderPlot({
    x <- current_obj(); req(x)
    p <- DimPlot(
      x,
      reduction = input$reduction,
      group.by  = if (isTRUE(input$setIdents)) NULL else input$groupby,
      label     = FALSE,
      pt.size   = 0.45,
      raster    = FALSE
    ) +
      ggtitle("UMAP (cosine)") +
      theme_classic(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", hjust = 0, margin = margin(b = 6)),
        axis.ticks      = element_blank(),
        legend.position = "bottom",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 10),
        legend.key.size = unit(8, "pt"),
        axis.line       = element_line(linewidth = 0.4)
      )
    p <- LabelClusters(plot = p,
                       id = if (isTRUE(input$setIdents)) "ident" else input$groupby,
                       repel = TRUE, size = 4)
    print(p)
  })
  
  # --- FeaturePlot ---
  output$p_feat <- renderPlot({
    x <- current_obj(); req(x)
    g <- genes(); req(length(g) >= 1)
    FeaturePlot(x, features = g, reduction = input$reduction, order = TRUE)
  })
  
  # --- Violin ---
  output$p_vln <- renderPlot({
    x <- current_obj(); req(x)
    g <- genes(); req(length(g) >= 1)
    VlnPlot(x, features = g, group.by = if (isTRUE(input$setIdents)) NULL else input$groupby) +
      theme(legend.position = "none")
  })
  
  # --- DotPlot ---
  output$p_dot <- renderPlot({
    x <- current_obj(); req(x)
    g <- genes(); req(length(g) >= 1)
    DotPlot(x, features = g, group.by = if (isTRUE(input$setIdents)) NULL else input$groupby) + RotatedAxis()
  })
  
  # --- マーカー（現在の Idents 基準） ---
  markersDf <- eventReactive(input$runMarkers, {
    x <- current_obj(); req(x)
    withProgress(message = "FindAllMarkers 実行中...", value = 0.1, {
      FindAllMarkers(x, only.pos = isTRUE(input$onlypos))
    })
  })
  
  output$tbl_markers <- renderDT({
    req(markersDf())
    datatable(markersDf(), options = list(pageLength = 10), rownames = FALSE)
  })
}

# ================= Run =================
shinyApp(ui, server)
