shinyPanelBatchcorrect <- fluidPage(
  tags$div(
    class = "container",
    h1("Batch Effect Diagnostics & Adjustment"),
    h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v06-tab04_Batch-Correction.html",
              "(help)", target = "_blank")),
    tabsetPanel(
        tabPanel(
          "Visualize Batch Effect",
          wellPanel(
            sidebarLayout(
              sidebarPanel(
                selectInput("batchAssay", "Select Assay:", currassays),
                tags$hr(),
                         
                selectInput("batchVarPlot", "Select Batch Annotation:", c("none", clusterChoice)),
                selectInput("conditionVarPlot", "Select Condition Annotation:", c("none", clusterChoice)),
                withBusyIndicatorUI(actionButton("visBatch", "Visualize batch effect"))#,
                #tags$hr(),
              ),
              mainPanel(
                wellPanel(
                  plotOutput("batchBoxplot", height = "600px")      
                )
              )
            )
          )
        ),
        
        tabPanel(
          "Run Batch Correction",
          wellPanel(
            sidebarLayout(
              sidebarPanel(
                selectInput("batchMethod", "Select Method:", c("ComBat", "ComBat-Seq")),
                tags$hr(),
                
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'ComBat'", "batchMethod"),
                  selectInput("combatAssay", "Select Assay:", currassays)
                ),
                
                selectInput("combatBatchVar", "Select Batch Variable:", clusterChoice),
                selectInput("combatConditionVar", "Select Primary Biological Variable:", clusterChoice),
                selectInput("combatCovariates", "Select Additional Covariates:", clusterChoice, multiple = TRUE),
                
                # original combat
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'ComBat'", "batchMethod"),
                  radioButtons("combatParametric", "Adjustments:", c("Parametric",
                                                                     "Non-parametric"),
                               selected = "Parametric"),
                  checkboxInput("combatMeanOnly", "Correct mean of the batch effect only",
                                value = FALSE),
                  checkboxInput("combatRef", "Run reference batch combat:",
                                value = FALSE),
                  uiOutput("selectCombatSeqBatchUI"),
                  textInput("combatSaveAssay", "Assay Name to Use:", value = "combat")
                ),
                
                # combat-seq
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'ComBat-Seq'", "batchMethod"),
                  checkboxInput("combatseqShrink", "Use Empirical Bayes Estimation:", value=FALSE),
                  conditionalPanel(
                    condition = sprintf("input['%s']", "combatseqShrink"),
                    textInput("combatseqGeneSubsetN", "Enter Number of Genes for Empirical Bayes:", "100")
                  ),
                  textInput("combatseqSaveAssay", "Assay Name to Use:", value = "combatseq")
                ),
                
                withBusyIndicatorUI(actionButton("combatRun", "Run batch correction"))    
              ),
              mainPanel(
                wellPanel(
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'ComBat'", "batchMethod"),
                    uiOutput("combatStatus"),
                    plotOutput("combatBoxplot", height = "600px")     
                  ),
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'ComBat-Seq'", "batchMethod"),
                    uiOutput("combatseqStatus"),
                    plotOutput("combatseqBoxplot", height = "600px")   
                  )
                )
              )
            )
          )
        )
      )
  )
)
