# ════════════════════════════════════════════════════════════════════════
# genomic-insights-apps/[nombre-analisis]/app.R
#
# Template de Shiny app para nuevos posts del blog.
# Copia esta carpeta, renómbrala y edita las secciones marcadas con ←
# ════════════════════════════════════════════════════════════════════════

library(shiny)
library(ggplot2)
library(dplyr)
library(bslib)

# ── Paleta corporativa del blog (no cambiar) ──────────────────────────────
VERDE       <- "#1D9E75"
VERDE_LIGHT <- "#E1F5EE"
AZUL        <- "#378ADD"
AMARILLO    <- "#FAC775"
ROSA        <- "#D4537E"

tema_app <- theme_minimal(base_family = "serif") +
  theme(
    plot.title       = element_text(size=13, face="bold", color="#1a1a1a"),
    plot.subtitle    = element_text(size=10, color="#666", margin=margin(b=8)),
    axis.text        = element_text(size=9, color="#444"),
    axis.title       = element_text(size=10, color="#333"),
    panel.grid.major = element_line(color="#ebebeb"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size=9),
    plot.background  = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="transparent", color=NA)
  )

# ════════════════════════════════════════════════════════════════════════
# ← DATOS: sustituye por la carga de tus archivos reales
# ════════════════════════════════════════════════════════════════════════
# datos <- readr::read_tsv("data/mis_datos.tsv")

datos_ejemplo <- data.frame(
  grupo = rep(c("A","B","C"), each = 20),
  valor = c(rnorm(20,3,0.5), rnorm(20,2,0.4), rnorm(20,4,0.6))
)

# ════════════════════════════════════════════════════════════════════════
# UI  ← añade tantas nav_panel() como figuras interactivas necesites
# ════════════════════════════════════════════════════════════════════════
ui <- page_navbar(
  title = "← Nombre del análisis",           # ← cambia
  theme = bs_theme(
    bootswatch   = "flatly",
    primary      = VERDE,
    base_font    = font_google("Source Serif 4"),
    heading_font = font_google("Playfair Display")
  ),
  bg = "#ffffff",

  nav_panel("Figura 1",                        # ← cambia el nombre
    layout_sidebar(
      sidebar = sidebar(
        title = "Controles",
        # ← añade aquí tus inputs Shiny:
        # selectInput("mi_input", "Etiqueta:", choices = c(...))
        # checkboxGroupInput(...)
        # sliderInput(...)
        selectInput("grupo", "Grupo:",
          choices = c("Todos","A","B","C"), selected = "Todos")
      ),
      card(
        full_screen = TRUE,
        card_header("← Título de la figura"),
        plotOutput("p1", height = "460px"),
        card_footer(
          class = "text-muted small fst-italic",
          "← Nota metodológica de la figura."
        )
      )
    )
  )

  # ← Duplica nav_panel() para añadir más figuras
)

# ════════════════════════════════════════════════════════════════════════
# SERVER
# ════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  output$p1 <- renderPlot({
    dat <- if (input$grupo == "Todos") datos_ejemplo
           else datos_ejemplo %>% filter(grupo == input$grupo)

    ggplot(dat, aes(x=grupo, y=valor, fill=grupo)) +
      geom_boxplot(alpha=0.4) +
      geom_jitter(alpha=0.5, width=0.15) +
      scale_fill_manual(values=c("A"=VERDE,"B"=AZUL,"C"=AMARILLO)) +
      labs(title="← Título", subtitle="← Subtítulo",
           x=NULL, y="← Eje Y") +
      tema_app + theme(legend.position="none")
  })

  # ← Añade más output$p2, output$p3... para cada figura adicional
}

shinyApp(ui, server)
