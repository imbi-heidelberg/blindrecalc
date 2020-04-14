#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (!require("shiny"))
  install.packages("shiny")
if (!require("devtools"))
  install.packages("devtools")
if (!require("blindrecalc"))
  devtools::install_github("imbi-heidelberg/blindrecalc")

library(shiny)
library(blindrecalc)

ui <- fluidPage(

   fluidPage(
     titlePanel("Blinded Sample Size Recalculation"),

     sidebarLayout(
       sidebarPanel(
                   sliderInput("alpha",
                               "Maximal Type One Error Rate",
                               min = .01,
                               max = .1,
                               step = .005,
                               value = .025),
                   sliderInput("beta",
                               "Maximal Type Two Error Rate",
                               min = .05,
                               max = .4,
                               step = .01,
                               value = .2),
                   sliderInput("r",
                               "Allocation Ratio",
                               min = 1/5,
                               max = 5,
                               step = .1,
                               value = 1),
                   sliderInput("delta",
                               "Alternative Effect Size",
                               min = .05,
                               max = 1,
                               value = .5,
                               step = .01),
                   sliderInput("delta_NI",
                               "Non-Inferiority Margin",
                               min = 0,
                               max = 1,
                               value = 0,
                               step = .01),
                   sliderInput("n_max",
                               "Maximal Sample Size",
                               min = 10,
                               max = 1000,
                               value = 100,
                               step = 1),
                   sliderInput("n1",
                               "Sample Size of Internal Pilot Study",
                               min = 1,
                               max = 1000,
                               value = 10,
                               step = 1),
                   sliderInput("nuisance",
                               "Assumed Variance",
                               min = 0.5,
                               max = 100,
                               value = 1,
                               step = .1),
                   sliderInput("iters",
                               "Number of Iterations",
                               min = 100,
                               max = 100000,
                               value = 1000,
                               step = 10),
                   submitButton("Compute!")

                  ),

          mainPanel(
            plotOutput(outputId = "dist_plot")#,
            )

          )
      ),

   p("This shiny app was developed by:"),
   a(href="https://www.klinikum.uni-heidelberg.de/medizinische-biometrie/wir-
ueber-uns/wir-ueber-uns", img(src="Logo-IMBI_ENGL.jpg", height=120))
)




server <- function(input, output) {

  output$dist_plot <- renderPlot({

    withProgress(message = "Computing characteristics",{

    d <- blindrecalc::setupStudent(alpha = input$alpha, beta = input$beta,
                                   r = input$r, delta = input$delta,
                                   delta_NI = input$delta_NI, n_max = input$n_max)

    return(blindrecalc::n_dist(d, input$n1, input$nuisance, summary = TRUE,
                               plot = TRUE, iters = input$iters))
    })

  })

}

# Run the application
shinyApp(ui = ui, server = server)

