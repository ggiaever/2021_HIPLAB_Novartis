
library(htmltools)
library(shiny)
library(dplyr)
library(shinyWidgets)
library(DT)
library(d3heatmap)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(ggplot2)
library(ggrepel)
library(shinyjs)


### RESPONSE SIGNATURES

### Response ALL signatures of combined UBC NIBR datasets
source(file.path("modules/2021_April25_AllSignature_Network.R"))

### Comparing response signatures from NIBR and UBC 
source(file.path("modules/2021_April23_Signature45_Network.R"))

source(file.path("modules/2021_April25_GO_visNetwork.R"))


source(file.path("modules/2021_April25_functions.R"))

source(file.path("modules/2021_April3_functions.R"))
source(file.path("modules/2021_April25_GOenrichment_functions.R"))

fdat    <-  read.delim("2021_April23_data/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)

xgold   <-  as.matrix(read.delim("2021_April23_data/2020_Sept3_aug30xgold.txt",
              stringsAsFactors = F, header = T, check.names = F))

xsig    <-  as.matrix(read.delim("2021_April23_data/xkq5_slim.txt",
              stringsAsFactors = F, header = T, check.names = F))



noness  <-  fdat  %>% filter(essential == "noness")
ess     <-  fdat  %>% filter(essential == "ess")



phiphop <- readRDS("2021_April23_data/2021_april25_phiphop.RDS")
xhiphop <- readRDS("2021_April23_data/2021_april13_xhiphop.RDS")

xcoinh  <- readRDS("2021_April23_data/2021_april13_coinhibition.RDS")
xpval   <- readRDS("2021_April23_data/2021_april13_coinhibition_pval.RDS")
xcofit  <- readRDS("2021_April23_data/2021_mar9_cofit.RDS")
xpfit   <- readRDS("2021_April23_data/2021_mar9_cofit_pval.RDS")
dend    <- readRDS("2021_April23_data/2020_oct28_aug30_dendrogram.RDS")
dend2   <- readRDS("2021_April23_data/2021_mar30_dendrogram.RDS")
we  <-  which(rownames(xhiphop)%in% ess$sgd_gene)
wn  <-  which(rownames(xhiphop)%in% noness$sgd_gene)

xe  <-  xhiphop[we,phiphop$name]
xne <-  xhiphop[wn,phiphop$name]

rm(noness)
rm(ess)

############
############


header <- dashboardHeader(title =  span("Chemogenomic Profiling HIPLAB vs NIBR",style = "font-weight:bold; font-size:24px"),titleWidth = 750)

sidebar <- dashboardSidebar(
  width = 250,
  tags$style(HTML(".main-sidebar {
                    background-color: #000066 !important; color: #ffffff !important; }
                   .treeview-menu>li>a { color: #ffffff !important;}")),
  sidebarMenu(

    id = "tabs", selected = "HIPHOP fitness profiles",

    menuItem(
      
       text = span("Gene fitness by screen", style = "font-size: 20px"),
       tabName="genebyscreen", icon = icon("chart-line")),

    menuItem(
      text = span("HIPHOP", style = "font-size: 20px"),
      tabName="hiphop", icon = icon("bullseye"), selected = T),

    menuItem(
      text = span("GO enrichments", style = "font-size: 20px"),
      tabName = "goenrich",icon = icon("dna")),
    
    menuItem(
      text = span("Response signature", style = "font-size: 20px"),
      tabName = "allsig",
      icon = icon("file-signature"),
      menuSubItem("All signatures", tabName = "signature"),
      menuSubItem("Compare signatures", tabName = "compsig"),startExpanded = T
    ),
    
    menuItem(
      text = span("Cluster analysis", style = "font-size: 20px"),
      tabName="heat",icon = icon("database"))
    
      )
  )



body <- dashboardBody(
  useShinyjs(), 
  
  #makes pngs same size as window after rescaling
  tags$script("$(document).on('shiny:connected', function(event) {
  var myWidth = $(window).width();
  Shiny.onInputChange('shiny_width',myWidth)

  });"),
  
  tags$script("$(document).on('shiny:connected', function(event) {
  var myHeight = $(window).height();
  Shiny.onInputChange('shiny_height',myHeight)

  });"),
  
   tags$head(
   tags$link(rel = "stylesheet", type = "text/css", href = "chemogenomics.css")
  ),
 
  
  tabItems(
    
    tabItem(tabName = "genebyscreen",
            
            fluidRow(
              column(width = 3,
                     
                     box(title = "Enter gene of interest:",
                         
                         textInput(
                           inputId = "gene",
                           label = "",
                           value   = "TOR2"
                         ), status = "primary", solidHeader = T, width = "100%", height = 250)
              ),
              
              column(width = 9,
                     box(title = "Gene descriptor and function:",
                         dataTableOutput('geneinfo')
                         ,status = "primary", solidHeader = T, width = "100%",height = 250)
              ),
              
              
              column(width = 12,
                     box(title = "Fitness defect score across screens:",
                         tags$div(
                           HTML("<h4><b>Click and drag to select screens of interest; profiles can be found by clicking on corresponding row 
                         in Experimental detail below AND are also listed in the compound menu on the HIPHOP tab. Graph: cyan: p-value < 0.001; 
                         triangle: indicates gene was identified as drug target in HIP screen.</b></h4>")),
                         br(),
                         
                         plotOutput("genebydrug", 
                                    brush = brushOpts(id = "brushpts"), height = 800),status = "primary", solidHeader = T
                         ,width = "80%")
              )
            ),
            
            fluidRow(
              
                box(title = "Experimental detail:",
                  
                  dataTableOutput("tabpts"),   
                  
                  status = "primary", solidHeader = T,width = 12)
              
            ),
            
            fluidRow(
              
                box(title="Cofit genes; select row to bring up gene profile:",
                  dataTableOutput('coFit'),status = "primary", solidHeader = T,width = 12)
            ),
            
            fluidRow(
                box(title="Download cofitness:", align = "center",
                  downloadButton('downloadcofitness', 'Download cofitness'),status = "primary", solidHeader = T,width = 3)
              
            )
    ),
    
    tabItem("hiphop",
            
            fluidRow(
              column(width = 6,
                     box(title = "Select HIPHOP screen:",status = "primary", solidHeader = T,width = "100%",height = 175,
                         selectizeInput("cmpSERV","", width = "400px",
                                        choices = NULL,
                                        multiple = F
                         )
                     )
              ),
              column(width = 3,
                     box(title = "Reset compound menu:",align = "center",height = 175,
                         br(),
                         actionButton(
                           inputId = "resetMenu",
                           label = "cmpMenu"),
                         
                         status = "primary",solidHeader = T, width = "100%")
                     
              )
              
            ),
            
            fluidRow(
              box(title = "HIP gene target & compound information:",status = "primary",
                  tags$div(
                    HTML("<h5><b>Click datatable row to view response signature:</b></h5>")),
                  
                  solidHeader = T, width = 12,DT::dataTableOutput("targethip")
              )),
            
            
            ######### HIP FITNESS PROFILES #########            
            fluidRow(
              
              box(title = "HaploIsufficient chemogenomic Profile (HIP):",status = "primary",
                  solidHeader = T,plotOutput("efd", width = "100%",height = 700)),
              
              box(title = "HOmozygous chemogenomic Profile (HOP):",status = "primary",
                  solidHeader = T,plotOutput("nfd", width = "100%",height = 700)),
            ),
            
            fluidRow(
              
              box(title="HIP genes",
                  dataTableOutput("genesHIP"),
                  status = "primary", solidHeader = T,width = 6),
              
              box(title="HOP genes",
                  dataTableOutput("genesHOP"),
                  status = "primary", solidHeader = T,width = 6)
            ),
            
            fluidRow(
              
              box(title="download HIP results:",
                  column(width = 6,
                         
                         downloadButton("downloadhip", "HIP download")
                    ),
                  
                  column(width = 6,
                         
                         downloadButton("ExportessPlot", "HIP fitness plot")
                    ),
                  
                  status = "primary", solidHeader = T,width = 6),
              
              
              box(title="download HOP results:",
                  column(width = 6,
                         
                         downloadButton("downloadhop", "HOP download")
                    ),
                  
                  column(width = 6,
                         
                         downloadButton("ExportnonPlot", "HOP fitness plot")
                    ),
                  
                  status = "primary", solidHeader = T,width = 6)
              
            ),
            
            fluidRow(
              column(width = 12,
                     box(title="Coinhibitory screens; select rows to bring up coinhibitory fitness profiles:",
                         dataTableOutput("coInhib"),status = "primary", solidHeader = T,width = "100%")
                    )
            ),
            
            fluidRow(
              box(title="download coinhibition:", align = "center",
                  downloadButton('downloadcoinhib', 'coinhibition'),status = "primary", solidHeader = T,width = 2)
            )
    ),
    tabItem("goenrich",
            fluidRow(
              visNetworkModuleUI("visNetwork1")   
            )
    ),
    
    tabItem("signature",
            fluidRow(
              sigNetworkModuleUI("sigNetworkModule1")
            ),
            
            fluidRow(
              box(title = "HIPHOP profiles in response signature:",
                  dataTableOutput("screensResp"),
                  status = "primary", solidHeader = TRUE,width = 12)
            )     
    ), 
    
    tabItem("compsig",    
            fluidRow(
              CompareResponseModuleUI("CompareResponseModule1")
            )
    ),
    
    tabItem("heat",
            fluidRow(
              box(title="Hierarchical cluster analysis of the combined UBC and NIBR datasets. To identify robust clusters, 
              we first compute coinhibition from the combined chemogenomic dataset, the pairwise Pearson correlation 
              between all screens to generate the “coinhibitory” square matrix, representing the similarity between profiled 
              compounds. Profiles were then hierarchically clustered using (1 – the coinhibitory matrix) as the distance 
              metric and the Ward agglomeration method. Branch colors on the two identical dendrograms indicate the major 
              clusters and are discussed below. Drugs within each cluster are highly correlated as indicated by the color scale as well 
              as by dendrogram height and correspond to compounds with similar chemogenomic profiles 
              that suggest similar mechanism of action.",
                  
                  d3heatmapOutput("heat", width=1800, height = 1600),solidHeader = T,
                  status = "primary", width = 12, height = 1600)
            ),
            
            fluidRow(
              box(title="Inspection of the resulting dendrograms reveal that the data primarily clusters by mechanism of action 
                  and not by research institute. Zoomed in clustering dendrogram clusters colored by (A) research institute; NIBR and UBC in navy and lightblue, respectively, 
                  and (B) by mechanism of drug action.",
                  plotOutput("dendy2", width = 1500, height = 600),
                  plotOutput("dendy", width = 1500, height = 600),solidHeader = T,
                  status = "primary", width = 12, height =1300)
            )
    )
  )
)
 ###########################################################################################################
  ui<-dashboardPage(header, sidebar, body,setBackgroundColor(color = "#e6e6ff", shinydashboard = TRUE))

 ###########################################################################################################
 
  server <- function(input, output, session){
    
  
  
##########################                          ##########################
##########################    CALL MODULES          ########################## 
##########################                          ##########################

    
#############################START MODULE FOR GO ENRICHMENT
tabSEND = reactive({
  
  if(is.null(input$tabs)) updateTabItems(session, "tabs", selected = "hiphop")
  
  input$tabs
  
})
  
xinput = reactive({
    xhiphop
    
  })

observeEvent(xinput(),{
  
  updateSelectizeInput(session,'cmpSERV',label = "",
                         choices = phiphop$name, selected = "NIBR_rapamycin:800pM", server = T)
})

cmpSEND = reactiveValues(cmp = NULL)
  
returnedCMP = callModule(visNetworkModule,'visNetwork1',xinput = reactive(xhiphop), 
                    cmp = reactive(cmpSEND$cmp),
                    
                    tab = tabSEND)
            

cmpRETURN = reactive({
                  
                  returnedCMP()})


############################ IMPORTANTANT 
### updates the SERVER compound if its different from the MODULE compound
### AND the user is currently on the GOENRICH fitness tab
############################ 
observeEvent(cmpRETURN(),{
            
  if(input$tabs == "goenrich" & (input$cmpSERV != cmpRETURN())){
    
    updateSelectizeInput(session, 'cmpSERV', choices = cmpRETURN(), server = T)
      
  }
},ignoreInit = F, ignoreNULL = T)


############################ IMPORTANTANT TO AVOID RECURSION
### updates the SERVER paramenters on any change;
### SENDS the SERVER compound only if the change occured on the HOP fitness tab
############################
observeEvent({input$cmpSERV
              input$tabs}, {
    
    choices = phiphop$name
    
    if(input$tabs == "hiphop"){
      cmpSEND$cmp = input$cmpSERV
      
      
      }
},ignoreInit = F, ignoreNULL = F)
#############################END MODULE FOR GO ENRICHMENT    

######################################################################  
######################################################################                                 
###########    START COMPARE SIGNATURE CALL MODULES          ######### 
######################################################################
######################################################################
 returnedCOMP <- callModule(CompareResponseModule,'CompareResponseModule1',

            inputTab  = tabSEND)

             
 
 
 
cntCOMPARE = reactive({
                  
                  returnedCOMP[[2]]()

                  })
 
 tabCOMPARE = reactive({
                  
                  returnedCOMP[[1]]()

                  })

 

# ################ receive the resp signature from the sigMODULE page 
# ################ SERV side drug signature table
  observeEvent({cntRETURN()
                tabRETURN()},{
    req(tabRETURN())
    
    newtab <- switch(input$tabs,
                     "compsig" = "signature",
                     "signature" = "compsig"
    )
    delay(500,updateTabItems(session, "tabs", newtab))  
   
   
  },ignoreInit = F, ignoreNULL = T)
  

  ################ SERV side drug signature table
  observeEvent({cntCOMPARE()
                tabCOMPARE()}, {
    req(tabCOMPARE())
    
    newtab <- switch(input$tabs,
                     "compsig" = "signature",
                     "signature" = "compsig")
                     

    delay(500,updateTabItems(session, "tabs", newtab))  
   
  },ignoreInit = F, ignoreNULL = T)
  
  
  
  
######################################################################  
######################################################################                                 
###########    END COMPARE SIGNATURE CALL MODULES          ########### 
######################################################################
######################################################################                              
                              
#######################################################################                              
#######################################################################                              
####  START ALL SIGNATURE MODULE FOR RESPONSE SIGNATURE ENRICHMENT ####
#######################################################################
####################################################################### 
  respSEND = reactiveValues(resp = NULL)

  

  returnedSIG <- callModule(sigNetworkModule,'sigNetworkModule1',xRespDat = reactive(xsig),
             xRespInput = reactive(respSEND$resp),

             inputTab  = tabSEND,

             message = "No GO enrichment, try another signature")

 
  sigRETURN = reactive({
                  
                  returnedSIG[[1]]()

                  })
  
  tabRETURN = reactive({
                  
                  returnedSIG[[2]]()

                  })
  

  cntRETURN = reactive({
               
                  returnedSIG[[3]]()
  
                  })

output$screensResp = renderDataTable({
  w = which(phiphop$signature %in% sigRETURN())

  validate(
    need(length(w) > 0 ,message =
           "please enter a valid signature"))
  d = phiphop$name[w]

  df = data.frame(screen = d,site = phiphop$site[w],mechanism = phiphop$mechanism[w],target = phiphop$target_html[w],
                  PCID = phiphop$pcid_link[w],FDA = phiphop$FDA[w],
                  gold_standard = phiphop$gold_standard[w], drug = phiphop$drug[w], stringsAsFactors = F)
  df = df %>% dplyr::arrange(desc(gold_standard),drug, site)
  w = which(names(df) %in% c("gold_standard","drug"))
  df = df[,-w]
  
  
  opts = list(dom = 'Bfrtip',
                autoWidth = T,scrollX = T,columnDefs = list(
                  list(className = 'dt-left',targets = c(0,2)),
                  list(className = 'dt-center',targets = c(1,3:5)),
                  list(width = c('50px'),targets = c(1,3,5)),
                  list(width = c('400px'),targets = c(2)),
                  list(width = c('50px'),targets = c(4)),
                  list(width = c('300px'),targets = 0)
                ))
    
    df =  DT::datatable(df, options = opts,rownames = F, escape = F, selection = "single") %>% 
      formatStyle(c(1:6),fontWeight = 'bold', fontSize = '14px')
        })
  


# ################ SERV side drug signature table from allSignature MODULE
  observeEvent(input$screensResp_rows_selected, {

    w = which(phiphop$signature %in% sigRETURN())
    d = phiphop[w,]
    row = input$screensResp_rows_selected


  validate(
    need(length(w) > 0 ,message =
           "please enter a valid signature"))

    d = phiphop$name[w]

    df = data.frame(screen = d,signature = phiphop$signature[w],target = phiphop$target[w],stringsAsFactors = F)
    df = df %>% arrange(screen)
    choices = df$screen[row]
    cmpSEND$cmp = choices

    

    updateSelectizeInput(session,"cmpSERV",label = "",
                         choices = df$screen[row], server = T)
                         
    newtab <- switch(input$tabs,
                     "hiphop" = "signature",
                     "signature" = "hiphop"
    )
     delay(500,updateTabItems(session, "tabs", newtab))

  })
################# send the signature to the signature MODULE page
 observeEvent(input$targethip_rows_selected, {
    newtab <- switch(input$tabs,
                     "hiphop" = "signature",
                     "signature" = "hiphop"
    )

    delay(500,updateTabItems(session, "tabs", newtab))
    w = which(phiphop$name %in% input$cmpSERV)
    d = phiphop[w,]

    row = input$targethip_rows_selected

    respSEND$resp = phiphop$signature[w][row][1]

  })
 
#######################################################################                              
#######################################################################                              
######  END ALL SIGNATURE MODULE FOR RESPONSE SIGNATURE ENRICHMENT ####
#######################################################################
####################################################################### 

 
 
    
 output$targethip <- DT::renderDataTable({
    w = which(phiphop$name %in% input$cmpSERV)
    d = phiphop[w,c("name","target_html","yeast_target_html","signature","mechanism","FDA","PCID")]
    d$FDA = toupper(d$FDA)
    d$FDA = factor(d$FDA)
    SGD ="<a href=https://www.yeastgenome.org/>known target</a>"
    PCID = "<a href=https://pubchem.ncbi.nlm.nih.gov/>pcid</a>"
    names(d)[1] = "screen"
    names(d)[2] = "HIP target"
    names(d)[7] = PCID
    
    names(d)[3] = "known target"
    d = d[,c("screen","HIP target","known target","signature","<a href=https://pubchem.ncbi.nlm.nih.gov/>pcid</a>","FDA","mechanism")]
    datatable(d,selection = "single", colnames = c("screen","HIP target","known target","signature","<a href=https://pubchem.ncbi.nlm.nih.gov/>pcid</a>","FDA","mechanism") ,



              options = list(dom = 'Bfrtip',paging = F,target = "cell",searching = F,info = F,
                    autowidth = T, scrollX = TRUE,
                    ordering = F,columnDefs = list(
                  
                  list(className = 'dt-left',targets = c(0,3,6)),
                  list(className = 'dt-center',targets = c(1,2,4,5)),
                  list(width = c('50px'),targets = c(2,4)),
                  list(width = c('75px'),targets = c(1)),
                  list(width = c('175px'),targets = c(0)),
                  list(width = c('300px'),targets = c(3)),
                  list(width = c('30px'),targets = c(5)),
                  list(width = c('250px'),targets = c(6))
                    )
                ),

              escape = F, class = 'table-bordered stripe table-condensed',rownames = F)
  })

                          
##########################                          ##########################
##########################    START GENEBYSCREEN    ########################## 
##########################                          ##########################
 
  data_for_genebydrug <- reactive({
    
    gene = toupper(input$gene)
    
    xdat = xhiphop
    
    w = which(rownames(xdat) == gene)
    
    validate(
      need(length(w) == 1 ,message = 
          "please enter a valid gene"))
    
    if(length(w) > 0) {

      mx2 = meltDF(xdat,row = gene,df = phiphop)
      
      mx2
    }
    
    mx2
  })
  
  
  ##########################
  ########################## PLOT GENEBYDRUG
  ##########################
  output$genebydrug <- renderPlot({
    x4 = data_for_genebydrug()
    mx2 = x4
    wna = which(is.na(mx2$fitness_defect))
    if(length(wna) > 0) mx2 = mx2[-wna,]
      
      tit = toupper(input$gene)
      mx2$fitness_defect = round(mx2$fitness_defect,2)
      mx2$sig = 0
      
      wx = which(mx2$fitness_defect > 5)
      if(length(wx) > 100) wx = wx[1:100]
      if(length(wx) > 0) mx2$sig[wx] = 1
      
      wsig = which(mx2$sig == T)
      
    
        g  = ggplot(mx2,aes(x =  screen,y=fitness_defect,col = factor(sig))) + theme_bw() +
             geom_point(aes(col = factor(sig)),shape = mx2$shape,size = 4) 
            

        g1 = g + theme(legend.position="none") +
          theme(panel.grid.minor =   element_blank()) +
          theme(panel.grid.major = element_blank()) +
          theme(axis.ticks = element_blank(), axis.text.x =   element_blank())

        g1 = g1 + labs(y = "fitness defect score") + labs(x = "compound") +
          geom_hline(yintercept=median(mx2$fitness_defect),color = "black",
            linetype = "dashed",size =1) +
          ggtitle(tit)
       
        g2 = g1 + geom_text_repel(size = 5,
          data = subset(mx2, sig == TRUE),
          aes(x =  screen,y=fitness_defect, label = screen),
          point.padding = 0.25,
          segment.alpha = 0.2, col = "black"
        )
    
        g2 = g2 + theme(axis.text=element_text(size=16),
          axis.title=element_text(size=24,face="bold"),
          plot.title = element_text(size = 24,face = "bold"))
  
        g2
      })
   
  ##########################
  ########################## GENEINFO GENEBYDRUG
  ##########################


    output$geneinfo  <- renderDataTable({
      w = which(fdat$sgd_gene %in% rownames(xhiphop))

      fd = fdat[w,] %>% distinct(sgd_orf,sgd_gene,.keep_all = T)

      w1 = which(fd$sgd_gene %in% toupper(input$gene))
      
      validate(
        need(length(w1) == 1 ,message =
               "please enter a valid gene"))
      
      f = fd[w1,c("sgd_orf","sgd_gene","descriptor","GOID","term")]
      
      
      f$sgd_gene = paste0("<a href=https://www.yeastgenome.org/locus/",f$sgd_gene,">",f$sgd_gene,"</a>")
      names(f)[4] = "GO ID"
      names(f)[5] = "GO biological process"
      names(f)[1:2] = c("ORF","GENE")
     
      f
    },escape = F,options=list(dom = 'Bfrtip',paging = F,searching = F, info = F,scrollX = T,autoWidth = F,
         columnDefs = list(
           list(width = "400px", targets = c(-3)),
           list(width = "100px", targets = c(-1)),
           list(width = "20px", targets = -c(2,4,5))
           )
         
         ),rownames = F)
    
    
##########################                          ##########################
##########################    START BRUSHPTS        ########################## 
##########################                          ##########################
 
  observeEvent(input$brushpts, {

    mx2 = data_for_genebydrug()
    
    pts = brushedPoints(mx2, input$brushpts, xvar = "screen", yvar = "fitness_defect")

    pts$fitness_defect = format(round( pts$fitness_defect, 2), nsmall = 2)
    
    g = grep("shape",names(pts))
    
    if(length(g) > 0) pts = pts[,-g]
    output$tabpts = DT::renderDataTable({


      pts


    },  options = list(dom = 'Bfrtip',paging=F,searching=F,info = F,
                    autowidth = T, scrollX = TRUE,
                    ordering=F,columnDefs = list(
                  list(className = 'dt-left',targets = -c(1:8)),
                  
                  list(width = c('100px'),targets = -c(2,5,8)),
                  
                  list(width = c('30px'),targets = -c(1,3,6,7)),
                  list(width = c('300px'),targets = -4)
                    )
                ),
              
              escape = F, class = 'table-bordered stripe table-condensed',rownames = F,
    
              selection = "single",server = F)
    

   choices = unique(pts$screen)

      updateSelectizeInput(session,"cmpSERV",label = "select condition",
                         choices = choices, selected = choices[1], server = T)


  },ignoreInit = T,ignoreNULL = F)

################# GENE BY SCREEN ROWS SELECTED #########  
  observeEvent(input$tabpts_rows_selected, {
      row = input$tabpts_rows_selected
      
      req(data_for_genebydrug())
      mx2 = data_for_genebydrug()
      
      nam = names(mx2)
     
      pts = brushedPoints(mx2, input$brushpts, xvar = "screen", yvar = "fitness_defect")
    
      choices = pts$screen[row]
      
      updateSelectizeInput(session,"cmpSERV",label = "select condition",
                         choices = choices, server = T)
      
      newtab <- switch(input$tabs,
                     "genebyscreen" = "hiphop",
                     "hiphop" = "genebyscreen"
                      )
      delay(500,updateTabItems(session, "tabs", newtab))
     
    })
 
##########################                          ##########################
##########################    START COFIT           ########################## 
##########################                          ##########################


outfit = reactive({
  w = which(colnames(xcofit) %in% toupper(input$gene))
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  p = xpfit[,w,drop=F]
  d = xcofit[,w,drop = F]
  
  d = d[order(d[,1],decreasing = T),,drop=F]
  
  df = data.frame(gene = rownames(d)[1:nrow(d)], cofitness = d[1:nrow(d),1],pvalue = p[rownames(d),1],stringsAsFactors = F)
  df = dfgeneAnno(df = df,gene = "gene",sgdlink = T)
  
  df$pvalue = formatC(df$pvalue,format = "e", digits = 2)
  
  df$cofitness = format(round( df$cofitness, 2), nsmall = 2)
  df = df[,c("ORF","GENE","cofitness","pvalue","description")]
  
  df
  
}
)                         
##########################                          
output$coFit = renderDataTable({
  fit = outfit()
  
  
},escape = F,selection = "single",rownames=F,options = list(
      autowidth = T,
      scrollX = TRUE,
      selection = "single",
      columnDefs = 
        list(
      list(width = '500px',targets = -1),
      list(width = '50px', targets = -c(5:2)),
      list(className = 'dt-center',targets = -c(4:2)),
      list(className = 'dt-left',targets = c(-5,-1))
          )))


 
  observeEvent(input$coFit_rows_selected, {
    row = input$coFit_rows_selected
    w = which(colnames(xcofit) %in% toupper(input$gene))
    validate(
      need(length(w) == 1 ,message = 
             "please enter a valid condition"))
    p = xpfit[,w,drop=F]
    d = xcofit[,w,drop = F]
    
    d = d[order(d[,1],decreasing = T),,drop=F]
    
    df = data.frame(gene = rownames(d)[1:500], cofit = d[1:500,1],pvalue = p[rownames(d)[1:500],1],stringsAsFactors = F)
    df$pvalue = formatC(df$pvalue,format = "e", digits = 2)
    
    df$cofit = format(round( df$cofit, 2), nsmall = 2)
    
    delay(500,updateTextInput(session, "gene", label = NULL, value = df$gene[row]))
    
    
  })

output$downloadcofitness <- downloadHandler(
  filename = function() {
    paste0("cofit:",toupper(input$gene), "_",Sys.Date(), ".txt")
    },
  content = function(file) {
    write.table(as.data.frame(outfit()), file, row.names = F,sep="\t",quote=F)
  }
)
##########################                          ##########################
##########################    END COFIT             ########################## 
   
##########################                            ##########################
##########################    END GENEBYSCREEN        ########################## 
##########################                            ##########################
  

  
##########################                          ##########################
##########################    START HIPHOP          ########################## 
##########################                          ##########################
##########################                          
observeEvent(input$resetMenu,{
    
    
    updateSelectizeInput(session,'cmpSERV',label = "",
                         choices = phiphop$name, selected = "NIBR_rapamycin:800pM", server = T)
    
    
  },ignoreInit = F, ignoreNULL = T)                       
##########################    GENE ANNOTATION OUTPUT                    
output$genesHIP = renderDataTable({
  req(input$cmpSERV)
  hip = geneAnno(mat = xe,fdat = fdat,sgdlink = T,cmp = input$cmpSERV)
  
  hip
},escape = F,
    options=list(pageLength=10, autoWidth = T,scrollX = TRUE,
                 columnDefs = list(list(width = "60px", targets = c(-2,-3)))),rownames = F)
  
  
  
  output$genesHOP = renderDataTable({
  req(input$cmpSERV)
  hop = geneAnno(mat = xne,fdat = fdat,sgdlink = T,cmp = input$cmpSERV)
   
  }, options=list(pageLength=10, autoWidth = T,scrollX = TRUE,
                 columnDefs = list(list(width = "60px", targets = c(-2,-3)))),rownames = F,escape = F)


##########################    FITNESS PLOT OUTPUT
 
  output$efd = renderPlot({
    validate(
      need(input$cmpSERV,message = 
          "please select condition"))
    cond1 = input$cmpSERV
    
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    wx = which(colnames(xe) %in% cond1)
    
    colnames(xe)[wx] = paste0("HIP|",colnames(xe)[wx])
    
    p10(xe,wx,pch = 17)
    
    
  })
  
  essPlot = function(){
    validate(
      need(input$cmpSERV,message = 
             "please select condition"))
    cond1 = input$cmpSERV
    
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    wx = which(colnames(xe) %in% cond1)
    
    colnames(xe)[wx] = paste0("HIP|",colnames(xe)[wx])
    
    p10(xe,wx,pch = 17)
    
  }
  
  output$nfd = renderPlot({
    validate(
      need(input$cmpSERV,message = 
             "please select condition"))
    cond1 = input$cmpSERV
    
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    wx = which(colnames(xne) %in% cond1)
    
    colnames(xne)[wx] = paste0("HOP|",colnames(xe)[wx])
    
    p10(xne,wx,pch = 17)
    
    
  })
  

  nonPlot = function(){
    validate(
      need(input$cmpSERV,message = 
             "please select condition"))
    cond1 = input$cmpSERV
    
    w = which(phiphop$name %in% cond1)
    exp = phiphop$name[w]
    
    wx = which(colnames(xne) %in% cond1)
    
    colnames(xne)[wx] = paste0("HOP|",colnames(xne)[wx])
    
    p10(xne,wx,pch = 17)
    
    
  }
  
  output$ExportessPlot <- downloadHandler(
    # file name
    filename = function() {
      paste0("HOPplot:",input$cmp1,  "_",Sys.Date(), ".png")
    },
    # content
    content = function(file){
      # create plot
      png(file,
          width = (input$shiny_width/12)*6,
          height = 1000)
      essPlot()
      dev.off()
    }) 
  
  output$ExportnonPlot <- downloadHandler(
    # file name
    filename = function() {
      paste0("HOPplot:",input$cmp1,  "_",Sys.Date(), ".png")
    },
    content = function(file){
      # create plot
      png(file,
          width = (input$shiny_width/12)*6,
          height = 1000)
      
      nonPlot()
      dev.off()
    })    
  
##########################    END FITNESS PLOT OUTPUT  

##########################
##########################
########################## GENE ANNOTATION OUTPUT

outhop = reactive({ 
  req(input$cmpSERV)
  hop = geneAnno(mat = xne,fdat = fdat,sgdlink = F,cmp = input$cmpSERV)
  
  })

output$downloadhop <- downloadHandler(
  filename = function() {
    paste0("HOP:",input$cmpSERV, "_",Sys.Date(), ".txt")
  },
  content = function(file) {
    write.table(as.data.frame(outhop()), file, row.names = F,sep="\t",quote=F)
  }
)  

outhip = reactive({
  req(input$cmpSERV)
  hip = geneAnno(mat = xe,fdat = fdat,sgdlink = F,cmp = input$cmpSERV)
 })


output$downloadhip <- downloadHandler(
  filename = function() {
    paste0("HIP:",input$cmpSERV, Sys.Date(), ".txt")
  },
  content = function(file) {
    write.table(as.data.frame(outhip()), file, row.names = F,sep="\t",quote=F)
  }
)
########################## END GENE ANNOTATION OUTPUT

##########################    END HIPHOP   ########################## 
##########################                 ##########################


##########################                          ##########################
##########################    START COINHIIBITION   ########################## 
##########################                          ##########################


output$coInhib = renderDataTable({
  
  xcoinhib = outcoInhib()
  
},escape = F,options = list(dom = 'Bfrtip', pageLength = 25),rownames=F,server = FALSE, selection = "single")


observeEvent(input$coInhib_rows_selected, {
    row <- input$coInhib_rows_selected
    df = outcoInhib()
   
    delay(300,updateSelectizeInput(session,"cmpSERV",label = "",
                                   choices = df$screen[row], selected = df$screen[row][1], , server = T))
    
  })


outcoInhib = reactive({
  w = which(colnames(xcoinh) %in% input$cmpSERV)
  validate(
    need(length(w) == 1 ,message = 
           "please enter a valid condition"))
  p = xpval[,w,drop=F]
  d = xcoinh[,w,drop = F]
  
  d = d[order(d[,1],decreasing = T),,drop=F]
  
  m = match(rownames(d),rownames(p))
  
  df = data.frame(screen = rownames(d)[1:500], coinhibition = d[1:500,1],pvalue = p[m[1:500],1],stringsAsFactors = F)
  df$pvalue = formatC(df$pvalue,format = "e",digits=2)
  df$coinhibition = format(round( df$coinhibition, 2), nsmall = 2)
  m = match(df$screen,phiphop$name)

  df$compound = phiphop$drug[m]
  df$known_target = phiphop$yeast_target_html[m]
  df$target = phiphop$target_html[m]
  df$signature = phiphop$signature[m]
  
  df$pcid = phiphop$pcid_link[m]
  df
  
}
)

output$downloadcoinhib <- downloadHandler(
  filename = function() {
    paste0("coinhibition:",input$cmpSERV, "_",Sys.Date(), ".txt")
    },
  content = function(file) {
    write.table(as.data.frame(outcoInhib()), file, row.names = F,sep="\t",quote=F)
  }
)
##########################    END COINHIIBITION   ########################## 
##########################                        ##########################

#######################                          ##########################
####################### DENDROGRAM AND HEATMAP PAGE
  output$heat = renderD3heatmap({
    col_fun = colorRamp2(c(0, 0.5, 1), c("midnightblue", "white", "darkred"))
    col = colorRampPalette(c("midnightblue","white","darkred"))(100)
    d3heatmap(xgold,Rowv = dend,Colv=dend2,width = 1200,height = 1000,yaxis_font_size = 10, xaxis_height = 380, yaxis_width = 380, xaxis_font_size = 10,
              col = col_fun,show_grid = T)
  })

  output$dendy = renderPlot({
    par(mar = c(25,4.1,1.1,2))
    plot(dend)
    
  })
  output$dendy2 = renderPlot({
    par(mar = c(25,4.1,1.1,2))
    plot(dend2)
  })
}


shinyApp(ui, server)
