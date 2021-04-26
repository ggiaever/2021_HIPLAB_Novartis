library(shiny)
library(dplyr)
library(DT)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(shinyjs)

source(file.path("modules","2021_April25_GOenrichment_functions.R"))
source(file.path("modules" ,"2021_April25_functions.R"))

fdat    <-  read.delim("2021_April23_data/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)
noness  <-  fdat  %>% filter(essential == "noness")
ess     <-  fdat  %>% filter(essential == "ess")


visNetworkModuleUI <- function(id, label = "visNetwork") {
   ns <- NS(id)
  
 tagList(
    fluidPage(title = "GO enrichments",
              

      
      ###########
      
      fluidRow(
        column(width  = 6,
               box(title = "Select HIPHOP screen:",
                   br(),
                   selectizeInput(ns('cmpMOD'),
                                  '', choices = NULL, multiple = F        
                   ),
                   status = "primary", solidHeader = T,width="100%", height = 200)
                ),
        
        
      
       column(width = 6, 
        tabBox(
          
          title = tagList(shiny::icon("info-circle"), tags$b("GO enrichment")),
          
          tabPanel("GO network",div(style="font-size:14px","Given a query set of genes (e.g. genes with significant FD scores
                  in a HIPHOP profile), we used the hypergeometric test to obtain
                  a P-value estimating the significance with which the set is
                  enriched with genes annotated to a given GO category. FDR values were obtained
                  using the Benjamini and Hochberg method. Analysis was restricted to genesets 
                  composed of > 5 and < 300 genes."),
         
          tags$a(href="https://www.ncbi.nlm.nih.gov/pubmed/24723613", "Lee et al. Science 2014.")),
          
          
          tabPanel("Top genes",div(style="font-size:14px","The enrichment of a GO gene set is driven by the gene subset
                        also present in the query set derived from the genes passing the fitness score threshold. 
                        In some cases, a similar set of query set genes drives the enrichment of 
                        multiple closely related GO gene sets. Of those in the same lineage in the GO hierarchy, 
                        we identified the most significantly enriched category and pruned the others.")),
                   
                   width = NULL,height = 200)
        
                )
        ),
      
  fluidRow(
    
    column(width = 3,
           box(title = "Reset compound menu:", align = "center",
               br(),
               actionButton(ns("resetGO"),label = "cmpMenu"),
               
               solidHeader = TRUE, width = "100%",height = 150, status = "primary")#box
    ),
    
        column(width = 3,
             box(title= "Background genome:", align = "center",
                 prettyRadioButtons(ns("pool"), label = "Background genome:", 
                                    c("HIP" = "ess", "HOP" = "noness","HIPHOP" = "both"), outline = T, fill = F,
                                    status = "primary", shape = "square",bigger = T,
                                    selected = "both", inline = T),
                 solidHeader = T,status = "primary",width = "100%",height = 150)),
    
    
        column(width = 3,       
             box(title = "Set fitness score threshold:",#input score threshold = 3 ~ P < 0.001
                 sliderInput(ns("scorethreshMOD"),label = "", min = 0, max = 5,
                             value = 3, step = 1,ticks = TRUE),
                 status = "primary", solidHeader = T,width="100%", height = 150)
              ),#coln
        
        column(width = 3,           
                box(title = "Set FDR threshold:",
                  sliderInput(ns("fdrMOD"),label = "", min = 0, max = 0.5,
                              value = 0.2, step = 0.1, ticks = TRUE), 
                  status = "primary", solidHeader = T,width="100%", height = 150)
              )), 
            
    
  fluidRow(
    
            box(title = "GO enrichment network: select a node to view details",
                width = 8, status = "primary", solidHeader = TRUE, height = 800,
                visNetworkOutput(ns("network_proxy"),width = "100%",height = 652)),
            
            box(title = "GO term set enrichment details:",width = 4, 
                uiOutput(ns("goTable")),solidHeader = T,status = "primary",background = "navy", height = 250),
            
            
            box(title = "Top-contributing genes:",width = 4,
                br(),
                uiOutput(ns("leadingEdge")),solidHeader = T,status = "primary", height = 530)
        ),
  
  fluidRow(
    box(title = "GO enrichment table:", DT::dataTableOutput(ns("enrtable")),
        status = "primary", solidHeader = TRUE,width = 12)
        ),
  fluidRow(
    column(width = 3,
           box(title= "Download enrichments:", align = "center",
               br(),
               downloadButton(ns("enrichdownload"), "enrichments"),
               solidHeader = T,status = "primary", width="100%", height = 150)
          )
        )
  ))
}


 #########################################################################################################################################     

visNetworkModule = function(input,output,session, xinput,cmp,tabs, message = "No GO enrichment, try relaxing the FDR or scorethreshold"){
  
  cmpRECVD = reactive({cmp()})
  datRECVD = reactive({xinput()})
  tabRECVD = reactive({tabs()})
  
  datMOD = reactiveValues(dat = NULL)
 ########################################################################
 ################  START GO ENRICHMENT NETWORK CALCULATIONS  ############
 ########################################################################
  
  observeEvent({input$pool
                datRECVD()},
              {
    req(datRECVD())
    req(input$pool)   
    ns <- session$ns 
    
    
    x = datRECVD()
    
    
    we  <-  which(rownames(x)%in% ess$sgd_gene)
    wn  <-  which(rownames(x)%in% noness$sgd_gene)
    
    if(input$pool == "noness") { xinp= x[wn,]} else if  (input$pool == "ess") {xinp = x[we,]} else {xinp = x}
    
    
    datMOD$dat = xinp
    
    })
  ############################ 
  goenrich <- reactive({ 
    req(datRECVD())
    req(!is.null(datMOD$dat))
    ns <- session$ns 
    
    xgo = datMOD$dat
    
    w = which(colnames(xgo) %in% input$cmpMOD)
    
    validate(
      need(length(w) != 0, "Please select a compound")
              )
    
    thresh = input$scorethreshMOD
    
    w2 = which(xgo[,w] >= thresh)
    
    validate(need(length(w2)!=0, message = "No scores above threshold"))

    req(length(w2)!=0)
    
    df = compSCORE(mat = xgo,coln = colnames(xgo)[w],sig = thresh)
   
    curr_exp = "network"
    
    FDR = input$fdrMOD
    
    network = runGOENRICH(fdrThresh = FDR, curr_exp = "tst",score = df, bp_path = "2021_April5_GOBP_SGD.RDS",go_path = "2021_April5_GOID_GOBP_SGD.txt")
    
    validate(
      need(!is.null(network$enrichInfo), message = message)
    )
   
    enrichInfo = network$enrichInfo
    
    req(!(is.null(enrichInfo)))
    
    edgeMat = network$edgeMat
    
    return(network)
  })
    
net <- reactive({ 
  
    req(goenrich()$enrichInfo)
    
    enrich = goenrich()$enrichInfo
    
    edge = goenrich()$edgeMat
    
    vis = visSetup(enrichInfo = enrich,edgeMat = edge, fontsize = 20, fontface = "Courier")
    
    
    vis
    
    })
 ########################################################################
 ################  END GO ENRICHMENT NETWORK CALCULATIONS  ##############
 ########################################################################
  



 ########################################################################
 ################  START  GO NETWORK, PLOTS & DATATABLES  ###############
 ########################################################################


output$network_proxy <- renderVisNetwork({
    req(net()$nodes)
    vis = net()
    ns <- session$ns 
    n = net()$nodes
    req(n)
    
    w =  nrow(n)
    n <- n %>% arrange(term)
    
    names = n$id
    
    
    if(nrow(vis$edges)==0) {
      visNetwork(vis$nodes, width = "100%") %>% 
        visNodes(shadow=list(enabled = T,size = 25),borderWidth=1) %>%
        visOptions(
          
          highlightNearest = list(enabled = T, degree = 5, hover = T),
          
          
          nodesIdSelection = list(enabled = TRUE, values = names,
             style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),
          
          selectedBy = list(variable="FDR",
            style = 'width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
    
    
    else {visNetwork(vis$nodes, vis$edges, width = "100%") %>% 
        visNodes(shadow=list(enabled = T,size = 25)) %>%
        
        visOptions(
          
          highlightNearest = list(enabled = T, degree = 5, hover = T),
          
          
          nodesIdSelection = list(enabled = TRUE, values = names,
              style = 'width: 500px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;'),
          
          selectedBy = list(variable="FDR",
              style = 'width: 200px; height = 31px; font-size: 18px; color: #000066;border: 3px solid #4d88ff;')) %>%
        
        visIgraphLayout(type = "full") %>%
        
        visEvents(select = "function(nodes) {
                Shiny.onInputChange('current_node_id', nodes.nodes);
                ;}")
    }
  })
############################  
  output$goTable = renderUI({
    
    req(input$network_proxy_selected)
    
    ns <- session$ns 
    
    DT::dataTableOutput(ns("table1"))
    
  })
############################
############################
############################
############################
############################ 
############################ 

############################
############################
############################ 
  output$leadingEdge = renderUI({
    req(input$network_proxy_selected)
    ns <- session$ns
    plotOutput(ns("bar"), width = 300,height = hgt())
  })
############################ 
output$bar <- renderPlot({

  req(input$network_proxy_selected)

  vis = net()

  n = net()$nodes

  w = which(vis$nodes$id %in% c(input$network_proxy_selected))

  req(length(w) > 0)

  n = vis$nodes[w,]

  req(nrow(n)!=0)

  s6 = geneBARPLOT(n$overlapGenes)

  barplot(s6$score,names.arg = s6$gene,las=1,horiz=T,col="navy")

})
############################ 
hgt = reactive({

  req(input$network_proxy_selected)

  vis = net()

  n = net()$nodes

  w = which(vis$nodes$id %in% c(input$network_proxy_selected))

  req(length(w) > 0)

  n = vis$nodes[w,]

  validate(need(nrow(n)!=0, message = "click node for detail"))

  o = geneBARPLOT(n$overlapGenes)

  height = genebarHEIGHT(o)

  height = height*2


  height
})
############################  
############################ 
  output$table1 <- DT::renderDataTable({
    
    req(input$network_proxy_selected)
    
    vis = net()
    
    n = net()$nodes
    
    w = which(vis$nodes$id %in% c(input$network_proxy_selected))
    
    
    req(length(w) > 0)
    
    term = vis$nodes$label[w]
    
    nam = c("term","nGenes","geneSetFraction","FDR")
    
    m = match(nam,names(vis$nodes))
    
    n = vis$nodes[w,m]
    
    term = vis$nodes$label[w]
    
    nam = c("term","nGenes","geneSetFraction","FDR")
    
    m = match(nam,names(vis$nodes))
    
    names(vis$nodes)[m] = c("GO term","geneSet size","% of geneSet","FDR")
    
    n = vis$nodes[w,m]
    
    req(nrow(n)!=0)
    
    t = t(n[,2:4])
    
    datatable(t,width=220,caption = htmltools::tags$caption(term,
                                                            style = "caption-side: top; text-align: center; color:black;background:white;font-weight:bold;"),
              options=list(paging=F,scrollY=F,dom="t",scroller=F,searching=F,ordering=F,rowCallback = JS(
                "function(row, data) {",
                "for (i = 1; i < data.length; i++) {",
                "if (data[i]>1000 | data[i]<1){",
                "$('td:eq('+i+')', row).html(data[i].toExponential(1));",
                "}",
                "}",
                "}")),
              height = 400,colnames = "") %>%
      formatStyle( target = "row", color = "black",backgroundColor = "white",
                   columns = c(1,2),fontWeight = "bold")
  })  
############################  
  
  
############################ 
  observeEvent(input$network_proxy_selectedBy,
              
               
               {
                 req(input$network_proxy_selectedBy)
                 
                 n = net()$nodes
                 
                 w = which(n$FDR %in% as.numeric(input$network_proxy_selectedBy))
                 
                 id = n$id[w]
                 
                 ns <- session$ns 
                 visNetworkProxy(ns("network_proxy")) %>%
                   visSelectNodes(id = id)
                 
               }
  )

 ########################################################################
 ################  END  GO NETWORK, PLOTS & DATATABLES  #################
 ########################################################################




 #################################################
 #################################################
 ########## START GO ENCHICHMENT TABLE ###########
 #################################################
 ################################################# 
  enrid = reactive({
    
    req(net()$nodes)
    
    enrich = net()$nodes
    
    req(length(nrow(enrich)) > 0)
    
    row <- input$enrtable_rows_selected
    
    out = outenrich()
    
    id = out$id[row]
    
    id
    
    
  })

  
############################  
outenrich = reactive({
    
    req(input$cmpMOD)
    req(goenrich()$enrichInfo)
    
    enrich = goenrich()$enrichInfo
    
    w = which(names(enrich) %in% c("querySetFraction", "geneSetFraction" ,
                                   "foldEnrichment" , "P" , "FDR" ))
    enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")] = 
      format(round(enrich[,c("querySetFraction","geneSetFraction", "foldEnrichment")],2),nsmall = 1,scientific = F)
    
    enrich[,c("P", "FDR")] = 
      format(signif(enrich[,c("P", "FDR")],2),nsmall = 1,scientific = T)
    
    enrich
  }
  )  
 ############################  
  observeEvent(enrid(), {
    req(enrid())
    id = enrid()
    
    ns <- session$ns
    visNetworkProxy(ns("network_proxy")) %>%
      visSelectNodes(id = id)
  })
############################   

  
############################ 
output$enrtable = DT::renderDataTable({
    
    out = outenrich()
    w = which(names(out) %in% c("GOID","term","querySetFraction", "geneSetFraction" ,
                                "foldEnrichment" , "P" , "FDR","overlapGenes" ))
    
    out = out[,c("GOID","term","querySetFraction", "geneSetFraction" ,
                 "foldEnrichment" , "FDR","overlapGenes" )]
    
    
    opts = list(dom = 'Bfrtip', pageLength = 25,
    autoWidth = T,scrollX = T,columnDefs = list(
                  list(className = 'dt-left',targets = c(0,1,6)),
                  list(className = 'dt-center',targets = c(2:5)),
                  list(width = c('40px'),targets = c(2:4)),
                  list(width = c('300px'),targets = c(6)),
                  list(width = c('50px'),targets = c(0,5)),
                  list(width = c('400px'),targets = 1)
                    ))
    
    df =  DT::datatable(out, options = opts,rownames = F, escape = F, selection = "single") %>% 
      formatStyle(c(1:7),fontWeight = 'bold', fontSize = '14px')
 })
  
  
############################ 
output$enrichdownload <- downloadHandler(
    filename = function() {
      paste0("enrich:",input$cmpMOD,"_", Sys.Date(), ".txt")
    },
    content = function(file) {
      write.table(outenrich(), file, row.names = F,sep="\t")
    }
  )
 #################################################
 ########## END GO ENCHICHMENT TABLE #############
 #################################################
 
  
 ################################################
 ###### START UPDATES & ObserveEVENTS ###########
 ################################################ 
 ######POPULATE cmpMOD

observeEvent(datRECVD(),{
  
    choices = colnames(datRECVD())
    updateSelectizeInput(session,'cmpMOD',label = "",
                         choices = choices, selected =  "NIBR_rapamycin:800pM", server = T)
    
    updateSliderInput(session,'scorethreshMOD',label = "", min = 0, max = 5,
                      value = 3, step = 1)
    updateSliderInput(session,'fdrMOD',label = "", min = 0, max = 0.5,
                      value = 0.2, step = 0.1)
    updatePrettyRadioButtons(session,"pool", label = "",
                               choices = c("HIP" = "ess", "HOP" = "noness","HIPHOP" = "both"), 
                               selected = "both", inline = T)
    
  },ignoreInit = F, ignoreNULL = T)

#################################################################################################
observeEvent(input$resetGO,{
  
    choices = colnames(datRECVD())
    updateSelectizeInput(session,'cmpMOD',label = "",
                         choices = choices, selected =  "NIBR_rapamycin:800pM", server = T)
    
    updateSliderInput(session,'scorethreshMOD',label = "", min = 0, max = 5,
                      value = 3, step = 1)
    updateSliderInput(session,'fdrMOD',label = "", min = 0, max = 0.5,
                      value = 0.2, step = 0.1)
    updatePrettyRadioButtons(session,"pool", label = "",
                               choices = c("HIP" = "ess", "HOP" = "noness","HIPHOP" = "both"), 
                               selected = "both", inline = T)
    
  },ignoreInit = F, ignoreNULL = T)

############################ IMPORTANTANT TASK
### updates the MODULE compound if its different from the SERVER compound
### AND the user is currently on the HOP fitness tab
############################ 
observeEvent(cmpRECVD(),{
    req(tabRECVD())
    req(cmpRECVD())

   # choices
   cmpSERV = cmpRECVD()
   # module cmp
   drugMOD = input$cmpMOD
    
   tabSERV = tabRECVD()
    
   #notSYNC = drugMOD != cmpSERV
     
    ns <- session$ns
    req(tabSERV)
    req(length(tabSERV > 0))
    
    ### if compound received from the server then update the compound in the module
    
    if(tabSERV != "goenrich"){
      
      updateSelectizeInput(session, 'cmpMOD', choices = cmpSERV)

      updateSliderInput(session,'scorethreshMOD',label = "", min = 0, max = 5,
                        value = 3, step = 1)
      updateSliderInput(session,'fdrMOD',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.1)
      updatePrettyRadioButtons(session,"pool", label = "",
                               choices = c("HIP" = "ess", "HOP" = "noness","HIPHOP" = "both"), 
                               selected = "both", inline = T)
     
    }
    
    
  },ignoreInit = F, ignoreNULL = T)
############################ IMPORTANT
### updates the MODULE paramenters on any change;
### RETURNS the MODULE compound only the change occured on the GOENRICH fitness tab
############################ THIS DIRECTION IS WORKING
observeEvent(input$cmpMOD,{
    req(input$cmpMOD)
    req(tabRECVD())
    req(cmpRECVD())
    
    
    
    ns <- session$ns 
    
    updateSliderInput(session,'scorethreshMOD',label = "", min = 0, max = 5,
                        value = 3, step = 1)
    updateSliderInput(session,'fdrMOD',label = "", min = 0, max = 0.5,
                        value = 0.2, step = 0.1)
    updatePrettyRadioButtons(session,"pool", label = "",
                               choices = c("HIP" = "ess", "HOP" = "noness","HIPHOP" = "both"), 
                               selected = "both", inline = T)
    
    notSYNC = input$cmpMOD != cmpRECVD()
    tabMOD = tabRECVD() == "goenrich" 
    
    ### if cmp was changed on the MOD side, then send the compound back to the server
    
    if(tabMOD & notSYNC){
        
        returnMOD$cmp = input$cmpMOD
        
        
        }

  },ignoreInit = F, ignoreNULL = T)
   
 ################################################
 ###### END UPDATES & ObserveEVENTS #############
 ################################################ 
  
returnMOD = reactiveValues(cmp = NULL)
############################ 
  return(
    
      reactive({ returnMOD$cmp })
      
    )
  
}

