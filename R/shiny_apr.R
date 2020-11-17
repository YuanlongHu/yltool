#' Run the TCMshiny Application
#'
#' @param ... A series of options to be used inside the app.
#'
#' @export
#' @importFrom shiny shinyApp
#' @author Yuanlong Hu
#' @examples
#'    AssoAnalysis()

AssoAnalysis <- function(){

  shiny::shinyApp(ui = app_ui_apr,
          server = app_server_apr
          )
}




#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shinydashboard
#' @importFrom shiny fluidRow
#' @importFrom shiny column
#' @importFrom shiny fileInput
#' @importFrom shiny radioButtons
#' @importFrom shiny numericInput
#' @importFrom shiny sliderInput
#' @importFrom visNetwork visNetworkOutput
#' @noRd
app_ui_apr <- function(request) {

  sidebar <- shinydashboard::dashboardSidebar(

    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Description", icon = icon("syringe"), tabName = "Description",
                               badgeColor = "blue"),
      shinydashboard::menuItem("Apriori", tabName = "Apriori", icon = icon("prescription"))

    )
  )

  tabltem_Description <- shinydashboard::tabItem(tabName = "Description",
                                                 fluidRow(
                                                   column(width = 12,
                                                          fluidRow(
                                                            column(width = 3,
                                                                   shinydashboard::box(collapsible = TRUE, width = 12,
                                                                                       fileInput("file1", "Choose CSV File",
                                                                                                 multiple = TRUE,
                                                                                                 accept = c("text/csv",
                                                                                                            "text/comma-separated-values,text/plain",
                                                                                                            ".csv")),
                                                                                       radioButtons("Format", "Format",
                                                                                                    choices = c(single = "single",
                                                                                                                basket = "basket"
                                                                                                    ),
                                                                                                    selected = "single"
                                                                                       )
                                                                   )
                                                            ),
                                                            column(width = 6,

                                                                   shinydashboard::tabBox(
                                                                     title = "Count", id = "Tab2",
                                                                     height = "600px", width = 12,
                                                                     DT::dataTableOutput("count")

                                                                   )

                                                            ),

                                                            column(width = 3,
                                                                   shinydashboard::box(collapsible = TRUE,
                                                                                       width = 12,
                                                                                       shinydashboard::valueBoxOutput("NumBox", width = 6),
                                                                                       shinydashboard::valueBoxOutput("HerbBox", width = 6)
                                                                   )
                                                            )
                                                          )
                                                   )
                                                 )
  )





  tabItem_Apriori <- shinydashboard::tabItem(tabName = "Apriori",
                                             fluidRow(
                                               column(width = 6,
                                                      fluidRow(
                                                        column(width = 12,
                                                               fluidRow(
                                                                 column(width = 12,
                                                                        shinydashboard::tabBox(
                                                                          title = "Results", id = "Tab",
                                                                          height = "600px", width = 12,
                                                                          DT::dataTableOutput("contents")

                                                                        )
                                                                 )
                                                               ),
                                                               fluidRow(
                                                                 column(width = 2,
                                                                        shinydashboard::box(collapsible = TRUE,width = 12,
                                                                                            numericInput("minlen", "minlen:", 2),
                                                                                            numericInput("maxlen", "maxlen:", 3)
                                                                        )
                                                                 ),

                                                                 column(width = 2,
                                                                        shinydashboard::box(collapsible = TRUE,width = 12,
                                                                                            radioButtons("Size", "Size",
                                                                                                         choices = c(Support = "support",Confidence = "confidence",Lift = "lift"),
                                                                                                         selected = "support"
                                                                                            )
                                                                        )
                                                                 ),
                                                                 column(width = 8,
                                                                        shinydashboard::box(collapsible = TRUE,width = 12,
                                                                                            sliderInput("support", "Support:",
                                                                                                        min=0, max=1, value=0.3,step = 0.01
                                                                                            ),
                                                                                            sliderInput("confidence", "Confidence:",
                                                                                                        min=0, max=1, value=0.5,step = 0.01
                                                                                            ),
                                                                                            numericInput("lift", "Lift:", 1)
                                                                        )
                                                                 )
                                                               )
                                                        )
                                                      )
                                               ),
                                               column(width = 6, visNetwork::visNetworkOutput("network", height = "1000px")

                                               )

                                             )
  )



  body <- shinydashboard::dashboardBody(
    shinydashboard::tabItems(
      tabltem_Description,
      tabItem_Apriori
    )
  )

  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "TCMshiny"),
    sidebar,
    body
  )


}



#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shinydashboard
#' @importFrom shiny req
#' @importFrom shiny reactive
#' @importFrom magrittr %>%
#' @importFrom visNetwork renderVisNetwork
#' @importFrom utils read.csv
#' @noRd
#' @author Yuanlong Hu

app_server_apr <- function(input, output, session) {
  # List the first level callModules here
  data <- reactive({
    req(input$file1)
    if(input$Format=="single"){
      data1 <- read.csv(input$file1$datapath,sep = ",",header = T)
      data1 <- data1[!duplicated(paste0(data1[,1], data1[,2])),]

    }else{
      data1 <- read.csv(input$file1$datapath,sep = ",",header = F)
    }
    return(data1)
  })

  # Input: datarules ----
  dataInput <- reactive({
    data <- data()
    datarules <- arules_apriori(data = data, minlen = input$minlen, maxlen = input$maxlen)
    return(datarules)
  })

  # Input: rules select ----
  datasetInput <- reactive({
    da <- dataInput()
    da <- select_rule(data = da, support = input$support, confidence = input$confidence, lift = input$lift)
    return(da)
  })

  # Input: count ----
  countInput <- reactive({

    if(input$Format=="single"){
      s <- count_herb(data())
    }

    names(s) <- c('Item','Frequency',"Proportion(%)")
    return(s)
  })


  # output: rules ----
  output$contents <- DT::renderDataTable({
    dat <- datasetInput()
    colnames(dat) <- c("A", "B", "C", "support", "confidence", "lift", "count")
    dat$A <- paste(dat$A, dat$B, dat$C)
    dat <- dat[,-c(2,3)]
    names(dat) <- c("Rule", "Support", "Confidence", "Coverage", "Lift", "Count")
    dat$Confidence <- signif(dat$Confidence,3)
    dat$Lift <- signif(dat$Lift,3)
    return(dat)

  },
  extensions = 'Buttons',
  rownames = F,
  options = list(searchHighlight = TRUE,
                 dom = 'lBfrtip',
                 scrollX = TRUE,
                 fixedColumns = list(leftColumns = 2, rightColumns = 1),
                 buttons = c('copy', 'csv', 'excel','pdf', 'print'),
                 lengthMenu = c(10,20,50,100,200))
  )


  # output: count ----
  output$count <- DT::renderDataTable({
    countInput()
  },
  extensions = 'Buttons',
  rownames = F,
  options = list(#language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Chinese.json'),
    searchHighlight = TRUE,
    dom = 'lBfrtip',
    scrollX = TRUE,
    fixedColumns = list(leftColumns =2, rightColumns = 1),
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
    lengthMenu = c(10,20,50,100,200))

  )

  output$NumBox <- shinydashboard::renderValueBox({
    Num <- length(unique(data()[,1]))

    shinydashboard::valueBox(
      Num, "Prescription", icon = icon("list"),
      color = "purple"
    )
  })
  output$HerbBox <- shinydashboard::renderValueBox({
    if(is.null(data())){
      Num <- 0
    }else{
      Num <- length(unique(data()[,2]))
    }

    shinydashboard::valueBox(
      Num, "Herb", icon = icon("list"),
      color = "light-blue"
    )
  })

  output$network <- visNetwork::renderVisNetwork({
    plot_data <- datasetInput()
    network_plot(plot_data, size = input$Size)
  })
}




#' apriori method using arules package
#'
#' @name arules_apriori
#' @param data input data
#' @param minlen minlen
#' @param maxlen maxlen
#' @importFrom arules apriori
#' @importFrom methods as
#' @noRd
#' @author Yuanlong Hu

arules_apriori <- function(data, minlen, maxlen){

  data <- data[,c(1:2)]
  colnames(data) <- c("TID","item")
  data <- methods::as(split(data[,"item"], data[,"TID"]), "transactions")

  data <- arules::apriori(data, parameter = list(support = 0.05,
                                         confidence = 0.05,
                                         minlen = minlen,
                                         maxlen = maxlen)
  )

  return(data)
}


#' select rule
#'
#' @name select_rule
#' @param data input data
#' @param support min support
#' @param confidence min confidence
#' @param lift min lift
#' @importFrom arules inspect
#' @noRd
#' @author Yuanlong Hu


select_rule <- function(data, support, confidence, lift = 1){

  da <- data
  dat <- subset(da,da@quality$support >= support & da@quality$confidence >= confidence & da@quality$lift > lift)
  dat <- arules::inspect(dat)
  dat <- as.data.frame(dat)
  return(dat)

}

#' count
#'
#' @name count_herb
#' @param data input data
#' @noRd
#' @author Yuanlong Hu

count_herb <- function(data){
  data <- as.data.frame(data)
  data <- data[,c(1:2)]
  s <- as.data.frame(table(data[,2]))
  names(s) <- c("drug","Freq")

  s$A <- round(s$Freq/length(unique(data[,1])),4)
  s$A <- s$A*100
  s <- s[order(s$Freq,decreasing = T),]
  names(s) <- c('Item','Frequency',"Proportion(%)")
  return(s)
}

#' count
#'
#' @name count_herb
#' @param data input data
#' @param size size
#' @import visNetwork
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub
#' @noRd
#' @author Yuanlong Hu



network_plot <- function(data, size="support"){

  if(size == "support"){
    data <- data[,-c(2,5:8)]
  }
  if(size == "confidence"){
    data <- data[,-c(2,4,6:8)]
  }
  if(size == "lift"){
    data <- data[,-c(2,4:6,8)]
  }
  data$lhs <- stringr::str_sub(data$lhs, start =  2, end = -2)
  data$rhs <- stringr::str_sub(data$rhs, start =  2, end = -2)
  data$ID <- paste0("node", c(1:nrow(data)))

  data <- data[,c(1,2,4,3)]


  data0 <- NULL
  for (i in 1:nrow(data)) {
    x <- as.character(data[i,])

    x1 <- strsplit(x[1], split = ",")[[1]]
    x1 <- data.frame(from = x1,
                     to=rep(x[3], length(x1)))
    x2 <- data.frame(from = x[3],
                     to= x[2])
    data1 <- rbind(x1, x2)
    data0 <- rbind(data0, data1)
  }
  data0$arrows <- "to"

  nodes <- unique(c(data0$from, data0$to))

  if(size=="lift"){
    nodes1 <- data.frame(id=nodes[nodes %in% data$ID],
                         label=nodes[nodes %in% data$ID],
                         size=data[,4]*20)
  }

  nodes1 <- data.frame(id=nodes[nodes %in% data$ID],
                       label=nodes[nodes %in% data$ID],
                       size=data[,4]*60)


  nodes2 <- data.frame(id=nodes[!nodes %in% data$ID],
                       label=nodes[!nodes %in% data$ID],
                       size=40)

  nodes <- rbind(nodes1, nodes2)
  nodes$group <- ifelse(nodes$label %in% data$ID, "lightblue", "white")
  nodes$label <- ifelse(nodes$label %in% data$ID, "", nodes$label)


  network_plot <- visNetwork::visNetwork(nodes, data0) %>%
    visNetwork::visEdges(shadow = F,
                         arrows =list(to = list(enabled = TRUE, scaleFactor = 2)),
                         color = list(color = "lightblue", highlight = "red")) %>%
    visNetwork::visOptions(highlightNearest = TRUE, manipulation = TRUE) %>%
    visNetwork::visGroups(groupname = "lightblue", color = "lightblue") %>%
    visNetwork::visGroups(groupname = "white", color = "white", shape = "text", size=40) %>%
    visNetwork::visInteraction(navigationButtons = TRUE, keyboard = FALSE) %>%
    visNetwork::visIgraphLayout(layout = "layout_nicely") %>%
    visNetwork::visExport()

  return(network_plot)
}
