packages <- c("ANN2", "gtools", "BiocManager", "tidyr", "shiny")
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "http://cran.us.r-project.org")
  }
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
lapply(c("ANN2", "gtools", "tidyr", "Biostrings", "shiny"), library, character.only = TRUE)
options(warn = -1) 
source("helper_functions.R")

ui <- fluidPage(
  titlePanel("GAP"),
  p("*** We recommend running GAP from a terminal, this app is meant for minimal use only!"),
  tabsetPanel(
    tabPanel("Predict unknown-phenotype status",
             sidebarLayout(
               sidebarPanel(
                 textInput("pp_transcript_id", "Transcript ID", "ENSMUST00000059970"),
                 textInput("pp_path_status", "Path to Status File", "input/species.txt"),
                 textInput("pp_path_file", "Path to Sequence File", "input/sample-dataset.fa"),
                 textInput("pp_path_tree", "Path to Tree File (Optional)", "input/tree-features.csv"),
                 checkboxInput("pp_use_tree_features", "Use Tree Features", value = FALSE),
                 #textInput("pp_path_to_results", "Path to Results", "results/"),
                 actionButton("pp_predict", "Predict")
               ),
               mainPanel(
                 textOutput("pp_status")
               )
             )),
    tabPanel("Positional p-values",
             sidebarLayout(
               sidebarPanel(
                 textInput("pv_transcript_id", "Transcript ID", "ENSMUST00000059970"),
                 textInput("pv_path_status", "Path to Status File", "input/species.txt"),
                 textInput("pv_path_file", "Path to Sequence File", "input/sample-dataset.fa"),
                 textInput("pv_path_tree", "Path to Tree File (Optional)", "input/tree-features.csv"),
                 checkboxInput("pv_use_tree_features", "Use Tree Features", value = FALSE),
                 #textInput("pv_path_to_results", "Path to Results", "results/"),
                 actionButton("pv_predict", "Calculate")
               ),
               mainPanel(
                 textOutput("pv_status")
               )
             )),
    tabPanel("Predict genomic region(s)",
             sidebarLayout(
               sidebarPanel(
                 #textInput("pr_transcript_id", "Transcript ID", "ENSMUST00000059970"),
                 textInput("pr_path_status", "Path to Status File", "input/species.txt"),
                 textInput("pr_path_file", "Path to Sequence File", "input/sample-dataset.fa"),
                 textInput("pr_path_tree", "Path to Tree File (Optional)", "input/tree-features.csv"),
                 checkboxInput("pr_use_tree_features", "Use Tree Features", value = FALSE),
                 #textInput("pv_path_to_results", "Path to Results", "results/"),
                 actionButton("pr_predict", "Predict")
               ),
               mainPanel(
                 textOutput("pr_status")
               )
             ))
  )
)

server <- function(input, output, session) {
  observeEvent(input$pp_predict, {
    if (input$pp_transcript_id == "" || input$pp_path_status == "" || input$pp_path_file == "") {
      output$pp_status <- renderText("Please fill all required fields.")
      return(NULL)
    }
    
    path_source <- getwd()
    path_order <- 'input/order.txt'
    path_status <- paste0(path_source, "/", input$pp_path_status)
    path_file <- paste0(path_source, "/", input$pp_path_file)
    path_tree <- ifelse(input$pp_path_tree == "", 'input/tree-features.csv', input$pp_path_tree)
    path_to_results <- paste0(path_source, "/", "results/")
    use_tree_features <- input$pp_use_tree_features
    numcv <- 1
    num_species <- get_num_species(path_status)
    path_order <- paste0(path_source, "/", path_order)
    
    try({
      PredictSpecies(
        transcript_id = input$pp_transcript_id,
        path_status = path_status,
        path_file = path_file,
        path_tree = path_tree,
        path_to_results = path_to_results,
        use_tree_features = use_tree_features,
        num_species = num_species,
        numcv = numcv,
        mincv = 1,
        po = path_order
      )
      output$pp_status <- renderText("Phenotype Prediction complete.")
    }, silent = TRUE)
  })
  
  observeEvent(input$pv_predict, {
    if (input$pv_transcript_id == "" || input$pv_path_status == "" || input$pv_path_file == "") {
      output$pv_status <- renderText("Please fill all required fields.")
      return(NULL)
    }
    
    path_source <- getwd()
    path_order <- 'input/order.txt'
    path_status <- paste0(path_source, "/", input$pv_path_status)
    path_file <- paste0(path_source, "/", input$pv_path_file)
    path_tree <- ifelse(input$pv_path_tree == "", 'input/tree-features.csv', input$pv_path_tree)
    path_to_results <- paste0(path_source, "/", "results/")
    use_tree_features <- input$pv_use_tree_features
    numcv <- 1
    num_species <- get_num_species(path_status)
    path_order <- paste0(path_source, "/", path_order)
    
    try({
      PositionalPvals(
        transcript_id = input$pv_transcript_id,
        path_status = path_status,
        path_file = path_file,
        path_tree = path_tree,
        path_to_results = path_to_results,
        use_tree_features = use_tree_features,
        num_species = num_species,
        numcv = numcv,
        mincv = 1,
        po = path_order
      )
      output$pv_status <- renderText("Positional p-values computation complete.")
    }, silent = TRUE)
  })
  
  observeEvent(input$pr_predict, {
    if (input$pr_path_status == "" || input$pr_path_file == "") {
      output$pr_status <- renderText("Please fill all required fields.")
      return(NULL)
    }
    
    path_source <- getwd()
    path_order <- 'input/order.txt'
    path_status <- paste0(path_source, "/", input$pr_path_status)
    path_file <- paste0(path_source, "/", input$pr_path_file)
    path_tree <- ifelse(input$pr_path_tree == "", 'input/tree-features.csv', input$pr_path_tree)
    use_tree_features <- input$pr_use_tree_features
    numcv <- 1
    hidden_layers<-0
    path_to_results<-paste0(path_source,"/results/associated.txt")
    if (!file.exists(path_to_results)) {
    file.create(path_to_results)
    }
    num_species <- get_num_species(path_status)
    path_order <- paste0(path_source, "/", path_order)
    ids<-as.list(get_unique_transcript_ids(path_file))
    try({
      PredictRegions(
        gene_list = ids,
        path_status = path_status,
        path_file = path_file,
        path_tree = path_tree,
        path_order = path_order,
        path_to_results = path_to_results,
        use_tree_features=use_tree_features,
        hidden_layers = hidden_layers,
        num_species = num_species,
        mincv = 0
        )
      output$pr_status <- renderText("Genomic region prediction complete.")
    }, silent = TRUE)
  })
}

shinyApp(ui = ui, server = server)