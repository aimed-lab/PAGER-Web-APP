# ui.R

#source('userInterface/vulcanoPlot.R')
# textOutput, textInput, uiOutput(), tableOutput(),plotlyOutput()
includeCSS("www/main.css")
### ui design ###
sidebar <- dashboardSidebar(
  skin = "light",
  sidebarMenu(
    id = "tabs", 
    bs4SidebarHeader("Dataset"),
    menuItem(
     "Load", 
     tabName = "loadTab" # link to the tabItem
    ),
    menuItem(
     "Filter", 
     tabName = "filterTab" # link to the tabItem
    ),
    menuItem(
     "Enrich", 
     tabName = "enrichTab" # link to the tabItem
    ),
    menuItem(
      "Term analysis", 
      tabName = "TanalysisTab" # link to the tabItem
    ),
    bs4SidebarHeader("Help"),
    menuItem(
      "About Us",
      tabName = "contactTab" #, icon = icon("users")
    )
  )
)

### contact ###
contactTab = tabItem(
  tabName = 'contactTab',
  shiny::tabsetPanel(
    tabPanel(
      includeMarkdown('userInterface/contact.Rmd')
    )
  )
)
### load ###
loadTab = tabItem(
  tabName = 'loadTab',
  #class = " main",
  inlineCSS(
    list(
      '.lightpink' = "background-color: lightpink",
      '.red'   = "background-color: red", 
      "textarea" = 'text-align: center', 
      '#text3 ' = 'text-align: center',
      'li a' = 'margin-right: 10px;',
      '.content-wrapper .content' = 'margin-top: 100px;',
      '.tabbable' = 'margin-top: 20px;'
    )
  ),
  
  #style = "display:flex;height:20px;text-align:center;
  #       flex-direction: column; justify-content: center;",
  shiny::tabsetPanel(
    tabPanel(
      title = "Load demo",
      fluidRow(
        column(
          width = 3,
          selectInput(
            "diseaseSelect", label = "Select a data set for analysis",
            choices = c('Melanoma',
                        'UALCAN_BLCA',
                        'UALCAN_BRCA',
                        'UALCAN_COAD',
                        'UALCAN_ESCA',
                        'UALCAN_HNSC',
                        'UALCAN_KICH',
                        'UALCAN_KIRC',
                        'UALCAN_KIRP',
                        'UALCAN_LIHC',
                        'UALCAN_LUAD',
                        'UALCAN_LUSC',
                        'UALCAN_PRAD',
                        'UALCAN_READ',
                        'UALCAN_STAD',
                        'UALCAN_THCA',
                        'UALCAN_UCEC'
                        ), #'lung squamous','breast cancer (Triple negative)'
            multiple = FALSE
         )
        ),
        column(
          width = 2,
          div(
            style = "display:inline-block; float:right", 
            actionButton(inputId = 'loadbtn', label = 'Proceed', status = "success"))
        ),
      ),
      
      includeMarkdown('userInterface/load.Rmd'),
      #HTML(markdown::markdownToHTML('userInterface/load.Rmd'))
      HTML("ULCAN cancer type abbreviation list"),
      fluidRow(
        DTOutput('tableCancer_list')   
      ),
    ),
    tabPanel(
      title = "Upload DEG table ",
      fluidRow(
        column(
          width = 12,
          fileInput(
            inputId = "DEtable", 
            label = "Upload your differentially expressed gene table",
            multiple = FALSE,
            accept = c(".txt",'.csv'),
          ),
          prettyCheckbox(
            inputId = "fheader1",
            label = "File contains header",
            value = TRUE,
            status = "primary",
            shape = "curve"
          )
        ),
        selectInput(
          "symbolcolumn", label = "Select the column of gene symbol",
          choices = ('symbol'),
          multiple = FALSE
        ),
        selectInput(
          "log2FCcolumn", label = "Select the column of log2FC",
          choices = ('log2FC'),
          multiple = FALSE
        ),
        selectInput(
          "pvaluecolumn", label = "Select the column of p-value",
          choices = ('pvalue'),
          multiple = FALSE
        )       
      ),
      column(
        width = 2,
        div(
          style = "display:inline-block; float:right", 
          actionButton(inputId = 'uploadbtn', label = 'Proceed', status = "success"))
      )
    ),

  ),

)

### filter ###
filterTab = tabItem(
  tabName = 'filterTab',
  class = " main",
  fluidRow(
    class = "volcano",
    style = "display:flex;height:20px;text-align:center;
         flex-direction: column; justify-content: center;",
    column(2,
           div(style = "display:inline-block; float:left", 
               actionButton('filterbackbtn', label = 'Back', status = "success"))),          
    column(8,
           HTML("<h5>Volcano plot of deferentially expressed genes.</h5> <em></em>")
    ),
    column(2,
           div(style = "display:inline-block; float:right", 
               actionButton('filterbtn', label = 'Proceed', status = "success"))
    ),
  ),
  hr(),
  ### text display ###
  fluidRow(
    column(3,
           textOutput('diseaseSelect')
    )
  ),
  h2("Gene filtering in the volcano plot"),
  fluidRow(
    sliderInput('sliderl2FC',
                label = 'log2foldchange',
                value = 1.0,
                min = 0.0,
                max = 10,
                step = 0.1
    ),
    sliderInput('sliderp_cut',
                label = '-log2p.value',
                value = 4.3,
                min = 0.0,
                max = 100,
                step = 0.1
    )
  ),
  column(2,
         div(style = "display:inline-block; float:right"),
         actionButton(inputId = "apply_btn_filter", label = "apply"),
  ),
  em("Click so changes will take effect."),
  fluidRow(
    plotlyOutput('vulcanoPlot')
  ),
  
  h2("Gene table after filtering"),
  fluidRow(
    DTOutput('tableDEGs_updated')   
  ),
)

### enrich ###
enrichTab = tabItem(
  tabName = 'enrichTab',
  fluidRow(
    class = "enrich",
    style = "display:flex;height:20px;text-align:center;
         flex-direction: column; justify-content: center;",
    column(2,
           div(style = "display:inline-block; float:left", 
               actionButton('enrichbackbtn', label = 'Back', status = "success"))),          
    column(8,
           HTML("<h5>Enriched terms of the downstream analysis tools</h5> <em></em>")
    ),
    column(2,
           div(style = "display:inline-block; float:right", 
               actionButton('enrichbtn', label = 'Proceed', status = "success"))
    ),
    hr(),
  ),
  hr(),
  shiny::tabsetPanel(
    tabPanel(
      title = "PAGER ",
      HTML("<h5>Enriched genesets</h5> <em></em>"),
      flowLayout(
        selectInput(
          inputId = 'selectPAGERsource',
          label = 'source',
          choices = c('BioCarta','DSigDB','GAD','GeneSigDB',
                      'GOA','GTEx','GWAS Catalog','Isozyme',
                      'KEGG_2021_HUMAN','Microcosm Targets','mirTARbase','MSigDB',
                      'NCI-Nature Curated','NGS Catalog','Pfam','PharmGKB',
                      'PheWAS','Protein Lounge','Reactome_2021','Spike',
                      'TargetScan','WikiPathway_2021','HPA-normProtein','HPA-PathologyAtlas',
                      'HPA-CellAtlas','HPA-RNAcon','HPA-normRNA','HPA-GTEx',
                      'HPA-FANTOM5','HPA-TCGA','Cell','I2D',
                      'I2D','GeoMx Cancer Transcriptome Atlas','CellMarker'
          ),
          selected = c('KEGG_2021_HUMAN',"WikiPathway_2021","BioCarta","NCI-Nature Curated","Reactome_2021","Protein Lounge","Spike"),
          multiple = T
        ),
        numericInput('sliderPAGERovlp',
                    label = 'Overlap>=',
                    value = 1,
                    min = 0,
                    max = 100,
                    step = 1
        ),         
        sliderInput('sliderPAGERsimilarity',
                    label = 'similarity score>=',
                    value = 0,
                    min = 0.0,
                    max = 1,
                    step = 0.01
        ),    
        sliderInput('sliderPAGERnlogpval',
                           label = '-log2p.value',
                           value = 4.3,
                           min = 0.0,
                           max = 100,
                           step = 0.1
        ),       
      ),  
      column(2,
        div(style = "display:inline-block; float:right"),
        actionButton(inputId = "apply_btn", label = "apply"),
      ),
      em("Click so changes will take effect."),
      fluidRow(
        plotlyOutput("PAGERhist"),height = NULL,width = NULL,status = "auto"
      ),
      fluidRow(
        plotlyOutput("PAGERpie"),height = NULL,width = NULL,status = "auto"
      ),      
      fluidRow(
        DTOutput('tablePAGER')   
      ),      
      fluidRow(
        HTML("<h5>M-type (comembership) PAG-to-PAG relationship</h5> <em></em>"),
        DTOutput('mtype_table')   
      ),
      HTML("<h5>M-type (comembership) PAG-to-PAG network</h5> <em></em>"),
      fluidRow(
        plotlyOutput('mtype')   
      ),

    ),
    tabPanel(
      title = "enrichR ", 
      HTML("<h5>Enriched genesets</h5> <em></em>"),
      flowLayout(
        selectInput(
          inputId = 'selectenrichRsource',
          label = 'source',
          choices = listEnrichrDbs()$libraryName[order(listEnrichrDbs()$libraryName)],
          selected = c("KEGG_2021_Human", "WikiPathway_2021_Human",
                       "BioCarta_2016","NCI-Nature_2016","Reactome_2016",
                       "Panther_2016"
                       ),
          multiple = T
        ),
        numericInput('sliderenrichRovlp',
                     label = 'Overlap>=',
                     value = 1,
                     min = 0,
                     max = 100,
                     step = 1
        ),         
        sliderInput('sliderenrichRnlogpval',
                    label = '-log2p.value',
                    value = 4.3,
                    min = 0.0,
                    max = 100,
                    step = 0.1
        ),       
      ),  
      column(
        2,
        div(style = "display:inline-block; float:right"),
        actionButton(inputId = "enrichR_apply_btn", label = "apply"),
      ),
      em("Click so changes will take effect."),
      fluidRow(
        plotlyOutput("enrichRhist"),height = NULL,width = NULL,status = "auto"
      ),
      fluidRow(
        DTOutput('tableenrichR_GS')   
      ),
    ),
    tabPanel(
      title = "WebGestaltR", 
      HTML("<h5>Enriched genesets</h5> <em></em>"),
      flowLayout(
        selectInput(
          inputId = 'selectwebgestaltRsource',
          label = 'source',
          choices = WebGestaltR::listGeneSet()$name,
          selected = c("pathway_KEGG","pathway_Reactome","pathway_Wikipathway"
                       ,"pathway_Panther"
                       ),
          multiple = T
        ),
        numericInput('sliderwebgestaltRovlp',
                     label = 'Overlap>=',
                     value = 1,
                     min = 0,
                     max = 100,
                     step = 1
        ),         
        sliderInput('sliderwebgestaltRnlogpval',
                    label = '-log2p.value',
                    value = 4.3,
                    min = 0.0,
                    max = 100,
                    step = 0.1
        ),       
      ),  
      column(2,
             div(style = "display:inline-block; float:right"),
             actionButton(inputId = "webgestaltR_apply_btn", label = "apply"),
      ),
      em("Click so changes will take effect."),
      fluidRow(
        plotlyOutput("webgestaltRhist"),height = NULL,width = NULL,status = "auto"
      ),
      fluidRow(
        DTOutput('tablewebgestaltR')   
      )
    ),
  ),
)
### analysis ###
TanalysisTab = tabItem(
  tabName = 'TanalysisTab',
  fluidRow(
    class = "Tanalysis",
    style = "display:flex;height:20px;text-align:center;
         flex-direction: column; justify-content: center;",
    column(2,
           div(style = "display:inline-block; float:left", 
               actionButton('Tanalysisbackbtn', label = 'Back', status = "success"))),       
    column(8,
           HTML("<h5>Term analysis result</h5> <em></em>")
    ),
    hr(),
  ),
  hr(),
  fluidRow(
    column(width = 12,
           sliderInput(
             inputId = 'sliderSimilairtyScore',
             label = 'Similarity score>=',
             value = 0.8,
             min = 0.6,
             max = 1.0,
             step = 0.01
           )
    ),

  ),
  column(2,
         div(style = "display:inline-block; float:right"),
         actionButton(inputId = "apply_btn_similarity", label = "apply"),
  ),
  em("Click so changes will take effect."),
  fluidRow(
    h5("Venn diagram."),
    plotOutput("venn")
    ##tabPanel("tabDAVID", plotOutput("DAVIDhist"))
  ),
  
  fluidRow(
    h5("Unique terms detected in PAGER."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(

        DTOutput('PAGER_uniq')
      )
    )
  ),  
  fluidRow(
    h5("Unique terms detected in EnrichR."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        
        DTOutput('enrichR_uniq')
      )
    )
  ), 
  fluidRow(
    h5("Unique terms detected in WebGestaltR."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        
        DTOutput('webgestaltR_uniq')
      )
    )
  ), 
  fluidRow(
    h5("Similar terms detected between WebGestaltR and PAGER."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        #title = "W_vs_P_uniq",
        DTOutput('W_vs_P_uniq')
      )
    )
  ),   
  fluidRow(
    h5("Similar terms detected between WebGestaltR and EnrichR."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        #title = "W_vs_E_uniq",
        DTOutput('W_vs_E_uniq')
      )
    )
  ),  
  fluidRow(
    h5("Similar terms between PAGER and EnrichR."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        #title = "P_vs_E_uniq",
        DTOutput('P_vs_E_uniq')
      )
    )
  ),
  fluidRow(
    h5("Similar terms among the PAGER, EnrichR and WebGestaltR."),
  ),
  fluidRow(
    shiny::tabsetPanel(
      tabPanel(
        #title = "Consensus",
        DTOutput('Consensus')
      )
    )
  ),
  
)

### body ###
body <- dashboardBody(
  tabItems(
    loadTab,
    filterTab,
    enrichTab,
    TanalysisTab,
    contactTab
  #tabItem(
  #  tabName = 'histomTab',
  #  box(
  #      plotlyOutput("plotTSNE", width = "auto", height = "auto"),
  #      height = NULL,
  #      width = NULL,
  #      status = "auto"
  #  )
  #),
  ),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "css/style.css")
  ),
  useShinyjs(),
)
### end body ###


### header ###
title <- dashboardBrand(
  title = "Analysis",
  color = "primary",
  href = "",
  opacity = 1
)

header <- dashboardHeader(
                fixed = TRUE,
                border = TRUE,
                status = 'info',
                title = title,
                div(style = "margin-left:auto;margin-right:auto; text-align:left; color:black",
                    HTML('<h3><strong>PAGER Web App</strong></h3>'),
                    HTML('<h6>An online interactive genesets and network analysis application.</h6>')
                )
                ,
                sidebarIcon = shiny::icon("water")
)
### end header ###


ui <- 
  dashboardPage(
    freshTheme = create_theme(
      bs4dash_vars(
        navbar_light_color = "#040404"
      ),
      bs4dash_layout(
        main_bg = "#fffffc" 
      ),
      bs4dash_sidebar_light(
        bg = "#FFF",
        color = "#040404",
        hover_color = "#0C4767",
      ),
      bs4dash_status(
        primary = "#4281A4", danger = "#BF616A", success = '#2a9d8f', warning = '#F7B538', info = "#fffffc"
      ),
      #    success = '#57A773'
      #    2a9d8f
      bs4dash_color(
        blue = '#4281A4', 
        lime = '#EBEBEB',
        white = '#fffffc'
      )
    ),    
    header,
    sidebar,  
    body,
    footer = dashboardFooter(
      left = a(
        href = "",
        target = "_blank", " "
      ),
      
      right = a(
        href = "",
        target = "_blank", HTML("")#<img height=25 src='images/aimed.png'/>
      )
    ),
    fullscreen = TRUE, dark = NULL
  )
