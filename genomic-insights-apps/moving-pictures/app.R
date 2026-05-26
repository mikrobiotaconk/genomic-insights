# ════════════════════════════════════════════════════════════════════════
# Moving Pictures Microbiome Explorer — app.R CORREGIDO
# ════════════════════════════════════════════════════════════════════════

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(bslib)
library(phyloseq)
library(vegan)
library(picante)
library(dunn.test)

# ── Paletas ───────────────────────────────────────────────────────────────
pal_body_site <- c("gut"="#1D9E75","left palm"="#EF9F27",
                   "right palm"="#D4537E","tongue"="#378ADD")
pal_subject   <- c("subject-1"="#1D9E75","subject-2"="#378ADD")
pal_antibiotic<- c("Yes"="#D4537E","No"="#1D9E75")

get_pal <- function(var) switch(var,
                                body_site  = pal_body_site,
                                subject    = pal_subject,
                                antibiotic = pal_antibiotic)

get_legend_title <- function(var) switch(var,
                                         body_site  = "Sitio corporal",
                                         subject    = "Sujeto",
                                         antibiotic = "Antibióticos")

get_labels <- function(var) switch(var,
                                   body_site  = c("gut"="Gut","left palm"="Piel izq.",
                                                  "right palm"="Piel der.","tongue"="Boca"),
                                   subject    = c("subject-1"="Sujeto 1","subject-2"="Sujeto 2"),
                                   antibiotic = c("Yes"="Con antibióticos","No"="Sin antibióticos"))

tema_app <- theme_minimal(base_family = "Arimo") +
  theme(
    plot.title       = element_text(family="Momo Trust Display", size=13,
                                    face="bold", color="#1a1a1a"),
    plot.subtitle    = element_text(family="Arimo", size=10, color="#666",
                                    margin=margin(b=8)),
    axis.text        = element_text(family="Arimo", size=9,  color="#444"),
    axis.title       = element_text(family="Arimo", size=10, color="#333"),
    panel.grid.major = element_line(color="#ebebeb"),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(size=9),
    legend.text      = element_text(size=9),
    plot.background  = element_rect(fill="transparent", color=NA),
    panel.background = element_rect(fill="transparent", color=NA)
  )

# ════════════════════════════════════════════════════════════════════════
# DATOS
# ════════════════════════════════════════════════════════════════════════
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

physeq <- readRDS("data/physeq.rds")

# ── Metadatos ─────────────────────────────────────────────────────────────
meta <- sample_data(physeq) %>%
  data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  rename(body_site  = body.site,
         antibiotic = reported.antibiotic.usage,
         days       = days.since.experiment.start) %>%
  mutate(body_site = factor(body_site,
                            levels = c("gut","left palm","right palm","tongue")),
         subject   = factor(subject))

# ── Diversidad alfa ───────────────────────────────────────────────────────
otu <- as.matrix(otu_table(physeq))
if (taxa_are_rows(physeq)) otu <- t(otu)

alpha_df <- data.frame(
  sample_id         = rownames(otu),
  shannon           = vegan::diversity(otu, index="shannon"),
  simpson           = vegan::diversity(otu, index="simpson"),
  observed_features = vegan::specnumber(otu),
  faith_pd          = picante::pd(otu, phy_tree(physeq),
                                  include.root=FALSE)$PD
) %>% left_join(meta, by="sample_id")

# ── Diversidad beta ───────────────────────────────────────────────────────
ord_bc <- ordinate(physeq, method="PCoA", distance="bray")
ord_wu <- ordinate(physeq, method="PCoA", distance="wunifrac")
ord_uu <- ordinate(physeq, method="PCoA", distance="unifrac")

extraer_pcoa <- function(ord, metrica) {
  coords <- as.data.frame(ord$vectors[,1:2])
  colnames(coords) <- c("Ax1","Ax2")
  pct <- round(ord$values$Relative_eig[1:2]*100,1)
  coords %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(meta, by="sample_id") %>%
    mutate(metrica=metrica, pct1=pct[1], pct2=pct[2])
}

pcoa_bc <- extraer_pcoa(ord_bc, "bc")
pcoa_wu <- extraer_pcoa(ord_wu, "wu")
pcoa_uu <- extraer_pcoa(ord_uu, "uu")

# ── Taxonomía ─────────────────────────────────────────────────────────────
ps_phylum <- tax_glom(physeq, taxrank="Phylum")

tax_data <- psmelt(ps_phylum) %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100) %>%
  ungroup() %>%
  rename(sample_id = Sample) %>%
  mutate(body_site = factor(body.site,
                            levels=c("gut","left palm","right palm","tongue"))) %>%
  group_by(body_site, Phylum) %>%
  summarise(abundancia=mean(rel_abund), .groups="drop") %>%
  group_by(body_site) %>%
  mutate(phylum_clean = ifelse(
    Phylum %in% names(sort(tapply(abundancia,Phylum,mean),
                           decreasing=TRUE))[1:5], Phylum, "Otros")) %>%
  group_by(body_site, phylum_clean) %>%
  summarise(abundancia=sum(abundancia), .groups="drop") %>%
  rename(phylum=phylum_clean) %>%
  mutate(body_site=factor(body_site,
                          levels=c("gut","left palm","right palm","tongue")))

phyla_unicos <- unique(tax_data$phylum)
pal_phylum   <- setNames(
  c("#1D9E75","#378ADD","#FAC775","#F09595","#5DCAA5","#B4B2A9",
    "#9B59B6","#E67E22","#1ABC9C","#E74C3C")[seq_len(length(phyla_unicos))],
  phyla_unicos)

# ── Helpers ───────────────────────────────────────────────────────────────
sitio_lbl <- c("gut"="Gut","left palm"="Piel izq.",
               "right palm"="Piel der.","tongue"="Boca")

label_ind <- function(x) switch(x,
                                shannon="Shannon (H')", simpson="Simpson (1-D)",
                                observed_features="Riqueza observada (ASVs)", faith_pd="Faith's PD")

# ════════════════════════════════════════════════════════════════════════
# UI
# ════════════════════════════════════════════════════════════════════════
ui <- page_navbar(
  title = "Moving Pictures Explorer",
  header = tags$style(HTML("
    .navbar { background-color: #085041 !important; padding: 0 1rem; }
    .navbar-nav .nav-item .nav-link {
      color: white !important;
      font-size: 0.88rem;
      padding: 0.5rem 1rem !important;
      border-radius: 6px;
      white-space: nowrap;
      margin: 4px 2px;
      border: 1px solid rgba(255,255,255,0.4) !important;
      opacity: 1 !important;
    }
    .navbar-nav .nav-item .nav-link:hover {
      background-color: rgba(255,255,255,0.2) !important;
      border-color: white !important;
    }
    .navbar-nav .nav-item .nav-link.active {
      background-color: #1D9E75 !important;
      border-color: #1D9E75 !important;
      font-weight: 600;
    }
  ")),
  theme = bs_theme(
    bootswatch   = "flatly",
    primary      = "#1D9E75",
    base_font    = font_google("Arimo"),
    heading_font = font_google("Momo Trust Display")
  ),
  bg = "#ffffff",
  
  # ── Tab 1: Diversidad alfa ──────────────────────────────────────────
  nav_panel("Diversidad alfa",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controles",
                selectInput("ind_alpha", "Índice de diversidad:",
                            choices  = c("Shannon"="shannon","Simpson (1-D)"="simpson",
                                         "Riqueza observada"="observed_features","Faith PD"="faith_pd"),
                            selected = "shannon"),
                selectInput("color_alpha", "Colorear por:",
                            choices  = c("Sitio corporal"="body_site","Sujeto"="subject",
                                         "Uso de antibióticos"="antibiotic"),
                            selected = "body_site"),
                checkboxGroupInput("sit_alpha", "Sitios:",
                                   choices  = c("Gut"="gut","Boca"="tongue",
                                                "Piel izq."="left palm","Piel der."="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                checkboxGroupInput("suj_alpha", "Sujetos:",
                                   choices  = c("Sujeto 1"="subject-1","Sujeto 2"="subject-2"),
                                   selected = c("subject-1","subject-2")),
                radioButtons("geom_alpha", "Tipo de gráfico:",
                             choices  = c("Boxplot"="box","Violín"="violin","Puntos + media"="dots"),
                             selected = "box"),
                checkboxInput("puntos_alpha", "Mostrar puntos sobre caja", value=TRUE)
              ),
              # ← figura correcta
              card(
                full_screen = TRUE,
                card_header("Diversidad alfa por sitio corporal"),
                plotOutput("p_alpha", height="460px"),    # ← p_alpha
                card_footer(class="text-muted small fst-italic",
                            "Rarefacción: 1.080 lecturas.")
              ),
              card(
                card_header("📊 Interpretación estadística"),
                uiOutput("texto_alpha")
              )
            )
  ),
  
  # ── Tab 2: PCoA ────────────────────────────────────────────────────
  nav_panel("Diversidad beta (PCoA)",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controles",
                selectInput("metrica", "Métrica de distancia:",
                            choices  = c("Bray-Curtis"="bc","UniFrac ponderado"="wu",
                                         "UniFrac no ponderado"="uu"),
                            selected = "bc"),
                selectInput("color_beta", "Colorear por:",
                            choices  = c("Sitio corporal"="body_site","Sujeto"="subject",
                                         "Uso de antibióticos"="antibiotic"),
                            selected = "body_site"),
                checkboxGroupInput("sit_beta", "Sitios:",
                                   choices  = c("Gut"="gut","Boca"="tongue",
                                                "Piel izq."="left palm","Piel der."="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                checkboxGroupInput("suj_beta", "Sujetos:",
                                   choices  = c("Sujeto 1"="subject-1","Sujeto 2"="subject-2"),
                                   selected = c("subject-1","subject-2")),
                checkboxInput("elipses",   "Elipses de confianza 95%", value=TRUE),
                checkboxInput("etiquetas", "Etiquetar muestras",       value=FALSE)
              ),
              # ← figura correcta
              card(
                full_screen = TRUE,
                card_header("Ordenación PCoA"),
                plotOutput("p_pcoa", height="460px"),     # ← p_pcoa
                card_footer(class="text-muted small fst-italic",
                            "PERMANOVA: R² = 0.68, p = 0.001 (sitio corporal).")
              ),
              card(
                card_header("📊 Interpretación estadística"),
                uiOutput("texto_beta")
              )
            )
  ),
  
  # ── Tab 3: Temporal ────────────────────────────────────────────────
  nav_panel("Evolución temporal",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controles",
                selectInput("sit_temp", "Sitio corporal:",
                            choices  = c("Gut"="gut","Boca"="tongue",
                                         "Piel izq."="left palm","Piel der."="right palm"),
                            selected = "gut"),
                selectInput("ind_temp", "Índice:",
                            choices  = c("Shannon"="shannon","Simpson (1-D)"="simpson",
                                         "Riqueza observada"="observed_features","Faith PD"="faith_pd"),
                            selected = "shannon"),
                selectInput("color_temp", "Colorear por:",
                            choices  = c("Sujeto"="subject","Uso de antibióticos"="antibiotic"),
                            selected = "subject"),
                checkboxInput("loess",  "Línea de tendencia (loess)", value=TRUE),
                checkboxInput("ribbon", "Intervalo de confianza",     value=TRUE),
                hr(),
                helpText("Cada sujeto mantiene un 'set point' personal de diversidad.")
              ),
              # ← figura correcta
              card(
                full_screen = TRUE,
                card_header("Dinámica temporal de la diversidad"),
                plotOutput("p_temp", height="460px"),     # ← p_temp
                card_footer(class="text-muted small fst-italic",
                            "Eje X: días desde el inicio del experimento.")
              ),
              card(
                card_header("📊 Interpretación estadística"),
                uiOutput("texto_temp")
              )
            )
  ),
  
  # ── Tab 4: Composición taxonómica ──────────────────────────────────
  nav_panel("Composición taxonómica",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controles",
                checkboxGroupInput("sit_tax", "Sitios:",
                                   choices  = c("Gut"="gut","Boca"="tongue",
                                                "Piel izq."="left palm","Piel der."="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                radioButtons("orient_tax", "Orientación:",
                             choices = c("Vertical"="v","Horizontal"="h"), selected="v"),
                hr(),
                helpText("Top 5 phyla por sitio; el resto agrupado en 'Otros'.")
              ),
              # ← figura correcta
              card(
                full_screen = TRUE,
                card_header("Abundancia relativa de phyla por sitio corporal"),
                plotOutput("p_tax", height="460px"),      # ← p_tax
                card_footer(class="text-muted small fst-italic",
                            "Abundancia relativa media (ambos sujetos). Clasificación: Silva 138.")
              )
            )
  )
)

# ════════════════════════════════════════════════════════════════════════
# SERVER
# ════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # ── Alpha ─────────────────────────────────────────────────────────────
  output$p_alpha <- renderPlot({
    req(input$suj_alpha, input$sit_alpha)
    
    dat <- alpha_df %>%
      filter(subject %in% input$suj_alpha, body_site %in% input$sit_alpha) %>%
      mutate(sitio_f   = factor(sitio_lbl[as.character(body_site)],
                                levels=sitio_lbl[input$sit_alpha]),
             color_var = as.character(.data[[input$color_alpha]]))
    
    y <- sym(input$ind_alpha); ylb <- label_ind(input$ind_alpha)
    pal <- get_pal(input$color_alpha); lbl <- get_labels(input$color_alpha)
    ltit <- get_legend_title(input$color_alpha)
    
    p <- ggplot(dat, aes(x=sitio_f, y=!!y, color=color_var, fill=color_var))
    
    if (input$geom_alpha == "box") {
      p <- p + geom_boxplot(alpha=0.2, width=0.5, position=position_dodge(0.7),
                            outlier.shape=if(input$puntos_alpha) NA else 19)
      if (input$puntos_alpha)
        p <- p + geom_jitter(alpha=0.5, size=1.8,
                             position=position_jitterdodge(0.12, dodge.width=0.7))
    } else if (input$geom_alpha == "violin") {
      p <- p + geom_violin(alpha=0.2, position=position_dodge(0.7)) +
        geom_boxplot(alpha=0, width=0.1, position=position_dodge(0.7))
    } else {
      p <- p + geom_jitter(alpha=0.65, size=2.8,
                           position=position_jitterdodge(0.15, dodge.width=0.65)) +
        stat_summary(fun=mean, geom="crossbar", width=0.35,
                     position=position_dodge(0.65), linewidth=0.9)
    }
    
    p + scale_color_manual(values=pal, name=ltit, labels=lbl) +
      scale_fill_manual( values=pal, name=ltit, labels=lbl) +
      labs(title=paste("Diversidad alfa —", ylb),
           subtitle=paste("Coloreado por:", ltit), x=NULL, y=ylb) +
      tema_app + theme(legend.position="top")
  })
  
  # ── Texto alpha ───────────────────────────────────────────────────────
  output$texto_alpha <- renderUI({
    req(input$ind_alpha, input$color_alpha, input$suj_alpha, input$sit_alpha)
    
    dat <- alpha_df %>%
      filter(subject %in% input$suj_alpha, body_site %in% input$sit_alpha)
    req(nrow(dat) > 0, length(unique(dat$body_site)) > 1)
    
    y <- dat[[input$ind_alpha]]
    kw     <- kruskal.test(y ~ dat$body_site)
    kw_p   <- format.pval(kw$p.value, digits=3, eps=0.001)
    kw_chi <- round(kw$statistic, 2)
    kw_df  <- kw$parameter
    
    dunn_html <- NULL
    if (length(unique(dat$body_site)) > 2) {
      dunn <- dunn.test::dunn.test(y, dat$body_site, method="bonferroni",
                                   alpha=0.05, kw=FALSE, label=TRUE)
      dunn_df <- data.frame(
        Comparacion   = dunn$comparisons,
        Z             = round(dunn$Z, 3),
        p_ajustado    = format.pval(dunn$P.adjusted, digits=3, eps=0.001),
        significativo = ifelse(dunn$P.adjusted < 0.05, "✓", "—")
      )
      filas <- apply(dunn_df, 1, function(r) {
        color_sig <- if (r["significativo"]=="✓") "#085041" else "#aaa"
        tags$tr(
          tags$td(r["Comparacion"],   style="padding:4px 8px; font-size:0.82rem;"),
          tags$td(r["Z"],             style="padding:4px 8px; font-size:0.82rem; text-align:center;"),
          tags$td(r["p_ajustado"],    style="padding:4px 8px; font-size:0.82rem; text-align:center;"),
          tags$td(r["significativo"], style=paste0("padding:4px 8px; font-size:0.82rem;
                  text-align:center; font-weight:600; color:",color_sig))
        )
      })
      dunn_html <- tagList(
        tags$p(style="font-size:0.78rem; font-weight:500; color:#085041; margin:0.75rem 0 0.25rem;",
               "Comparaciones post-hoc (Dunn · corrección Bonferroni):"),
        tags$table(style="width:100%; border-collapse:collapse;",
                   tags$thead(tags$tr(style="border-bottom:1px solid #9FE1CB;",
                                      tags$th("Comparación", style="padding:4px 8px; font-size:0.78rem; text-align:left;"),
                                      tags$th("Z",           style="padding:4px 8px; font-size:0.78rem; text-align:center;"),
                                      tags$th("p ajustado",  style="padding:4px 8px; font-size:0.78rem; text-align:center;"),
                                      tags$th("Sig.",        style="padding:4px 8px; font-size:0.78rem; text-align:center;")
                   )),
                   tags$tbody(filas)
        )
      )
    }
    
    texto <- switch(input$color_alpha,
                    body_site = switch(input$ind_alpha,
                                       shannon           = tags$p("Shannon pondera riqueza y equitatividad. La cavidad oral lidera (H'~3.5), seguida del intestino (~3.0) y la piel (~2.0)."),
                                       simpson           = tags$p("Simpson da más peso a especies dominantes. Los patrones se mantienen pero con diferencias menos pronunciadas que Shannon."),
                                       observed_features = tags$p("La riqueza bruta de ASVs sigue el mismo gradiente pero es más sensible a la profundidad de secuenciación."),
                                       faith_pd          = tags$p("Faith's PD incorpora distancia evolutiva entre linajes — revela diversidad filogenética, no solo taxonómica."),
                                       tags$p("Selecciona un índice.")),
                    subject = switch(input$ind_alpha,
                                     shannon           = tags$p("La variación entre sujetos es menor que entre sitios. Cada individuo mantiene un 'set point' propio estable en el tiempo."),
                                     simpson           = tags$p("Con Simpson las diferencias entre sujetos son aún menos marcadas — dominancia de especies principales similar entre individuos."),
                                     observed_features = tags$p("La riqueza bruta varía más entre sujetos que los índices ponderados, siendo más sensible a diferencias en composición específica."),
                                     faith_pd          = tags$p("La diversidad filogenética entre sujetos es sorprendentemente similar pese a diferir en composición taxonómica."),
                                     tags$p("Selecciona un índice.")),
                    antibiotic = switch(input$ind_alpha,
                                        shannon           = tags$p("Los antibióticos reducen Shannon claramente, especialmente en intestino. La recuperación es visible en muestras posteriores."),
                                        simpson           = tags$p("Con Simpson el efecto es más visible: dominancia de pocas especies resistentes aumenta durante el tratamiento."),
                                        observed_features = tags$p("La riqueza cae de forma abrupta — se pierden ASVs enteros. Es el índice más sensible para detectar el impacto inicial."),
                                        faith_pd          = tags$p("Una caída de PD indica pérdida de linajes evolutivos enteros — más difícil de recuperar que la riqueza taxonómica."),
                                        tags$p("Selecciona un índice.")),
                    tags$p("Selecciona una variable."))
    
    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretación — ", get_legend_title(input$color_alpha),
                           " · ", label_ind(input$ind_alpha))),
        tags$div(style="margin-top:0.6rem; padding:0.5rem 0.75rem;
                      background:white; border-radius:6px;
                      font-size:0.82rem; color:#333;",
                 tags$span(style="font-weight:500;", "Kruskal-Wallis: "),
                 paste0("H(", kw_df, ") = ", kw_chi, ",  p = ", kw_p),
                 if      (kw$p.value < 0.001) tags$span(style="color:#1D9E75; font-weight:600; margin-left:6px;", "***")
                 else if (kw$p.value < 0.01)  tags$span(style="color:#1D9E75; font-weight:600; margin-left:6px;", "**")
                 else if (kw$p.value < 0.05)  tags$span(style="color:#1D9E75; font-weight:600; margin-left:6px;", "*")
                 else                          tags$span(style="color:#aaa; margin-left:6px;", "n.s.")
        ),
        dunn_html,
        tags$div(style="margin-top:0.75rem;", texto)
    )
  })
  
  # ── PCoA ──────────────────────────────────────────────────────────────
  output$p_pcoa <- renderPlot({
    req(input$sit_beta, input$suj_beta)
    
    dat <- switch(input$metrica, bc=pcoa_bc, wu=pcoa_wu, uu=pcoa_uu) %>%
      filter(body_site %in% input$sit_beta, subject %in% input$suj_beta) %>%
      mutate(sitio_f   = factor(sitio_lbl[as.character(body_site)],
                                levels=sitio_lbl[input$sit_beta]),
             color_var = as.character(.data[[input$color_beta]]))
    
    pct1 <- dat$pct1[1]; pct2 <- dat$pct2[1]
    nom  <- switch(input$metrica, bc="Bray-Curtis",
                   wu="UniFrac ponderado", uu="UniFrac no ponderado")
    pal  <- get_pal(input$color_beta)
    lbl  <- get_labels(input$color_beta)
    ltit <- get_legend_title(input$color_beta)
    
    p <- ggplot(dat, aes(x=Ax1, y=Ax2, color=color_var, fill=color_var)) +
      geom_point(size=4, alpha=0.85) +
      scale_color_manual(values=pal, name=ltit, labels=lbl) +
      scale_fill_manual( values=pal, name=ltit, labels=lbl, guide="none") +
      labs(title=paste0("PCoA (", nom, ")"),
           subtitle=paste("Coloreado por:", ltit),
           x=paste0("PCo1 (",pct1,"%)"), y=paste0("PCo2 (",pct2,"%)")) +
      tema_app + theme(legend.position="right")
    
    if (input$elipses)
      p <- p + stat_ellipse(aes(group=color_var), level=0.95,
                            linewidth=0.65, linetype="dashed", alpha=0.55)
    if (input$etiquetas)
      p <- p + geom_text(aes(label=color_var), size=2.5,
                         nudge_y=0.015, alpha=0.8)
    p
  })
  
  # ── Texto beta ────────────────────────────────────────────────────────
  output$texto_beta <- renderUI({
    req(input$metrica, input$color_beta)
    
    nom  <- switch(input$metrica, bc="Bray-Curtis",
                   wu="UniFrac ponderado", uu="UniFrac no ponderado")
    ltit <- get_legend_title(input$color_beta)
    
    texto <- switch(input$color_beta,
                    body_site = switch(input$metrica,
                                       bc = tags$p("Bray-Curtis separa claramente los cuatro sitios corporales. El sitio explica el 68% de la varianza total (PERMANOVA R²=0.68, p=0.001) — mucho más que el sujeto individual (R²=0.09)."),
                                       wu = tags$p("UniFrac ponderado incorpora abundancia relativa y distancia filogenética. La separación entre sitios se mantiene pero refleja diferencias evolutivas más profundas entre comunidades."),
                                       uu = tags$p("UniFrac no ponderado es más sensible a especies raras. Puede revelar diferencias entre sitios que Bray-Curtis no detecta por estar dominado por las especies más abundantes.")),
                    subject = switch(input$metrica,
                                     bc = tags$p("Con Bray-Curtis las muestras del mismo sujeto tienden a agruparse dentro de cada sitio corporal, pero la separación por sujeto es mucho menor que por sitio (R²=0.09 vs 0.68)."),
                                     wu = tags$p("UniFrac ponderado puede revelar diferencias entre sujetos más sutiles a nivel filogenético, aunque el patrón dominante sigue siendo el sitio corporal."),
                                     uu = tags$p("UniFrac no ponderado a veces amplifica diferencias entre sujetos al dar más peso a especies raras que pueden ser individuo-específicas.")),
                    antibiotic = switch(input$metrica,
                                        bc = tags$p("El uso de antibióticos puede desplazar las muestras en el espacio PCoA, especialmente en intestino. Las muestras bajo tratamiento tienden a alejarse del cluster habitual de su sitio."),
                                        wu = tags$p("UniFrac ponderado refleja bien el impacto de los antibióticos a nivel filogenético — la pérdida de linajes enteros desplaza las muestras más que la simple reducción de abundancia."),
                                        uu = tags$p("UniFrac no ponderado es especialmente sensible al efecto de antibióticos sobre especies raras, que pueden desaparecer completamente durante el tratamiento.")))
    
    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretación — ", ltit, " · ", nom)),
        tags$div(style="margin-top:0.75rem;", texto)
    )
  })
  
  # ── Temporal ──────────────────────────────────────────────────────────
  output$p_temp <- renderPlot({
    dat <- alpha_df %>%
      filter(body_site == input$sit_temp) %>%
      mutate(color_var = as.character(.data[[input$color_temp]]))
    
    y <- sym(input$ind_temp); ylb <- label_ind(input$ind_temp)
    nom  <- sitio_lbl[input$sit_temp]
    pal  <- get_pal(input$color_temp)
    lbl  <- get_labels(input$color_temp)
    ltit <- get_legend_title(input$color_temp)
    
    p <- ggplot(dat, aes(x=days, y=!!y, color=color_var, group=color_var)) +
      geom_line(linewidth=0.9, alpha=0.75) +
      geom_point(size=2.5, alpha=0.85) +
      scale_color_manual(values=pal, name=ltit, labels=lbl) +
      labs(title=paste0("Variación temporal — ", nom),
           subtitle=paste0(ylb, " · coloreado por: ", ltit),
           x="Días desde el inicio del experimento", y=ylb) +
      tema_app + theme(legend.position="top")
    
    if (input$loess)
      p <- p + geom_smooth(aes(group=color_var), method="loess",
                           se=input$ribbon, linewidth=0.55,
                           alpha=0.12, linetype="dashed")
    p
  })
  
  # ── Texto temporal ────────────────────────────────────────────────────
  output$texto_temp <- renderUI({
    req(input$sit_temp, input$ind_temp, input$color_temp)
    
    nom  <- sitio_lbl[input$sit_temp]
    ylb  <- label_ind(input$ind_temp)
    ltit <- get_legend_title(input$color_temp)
    
    texto <- switch(input$color_temp,
                    subject = switch(input$ind_temp,
                                     shannon           = tags$p("Cada sujeto mantiene un 'set point' propio de Shannon relativamente estable. Las fluctuaciones son moderadas y el sistema tiende a recuperar su nivel basal tras perturbaciones."),
                                     simpson           = tags$p("Con Simpson la estabilidad individual es aún más marcada — la dominancia de las especies principales es muy consistente en el tiempo dentro de cada individuo."),
                                     observed_features = tags$p("La riqueza observada fluctúa más que los índices ponderados, siendo más sensible a eventos puntuales como cambios de dieta o infecciones leves."),
                                     faith_pd          = tags$p("La diversidad filogenética muestra la mayor estabilidad temporal — los linajes evolutivos presentes cambian poco aunque cambien las especies individuales.")),
                    antibiotic = switch(input$ind_temp,
                                        shannon           = tags$p("El tratamiento antibiótico provoca caídas puntuales de Shannon visibles en el eje temporal. La velocidad de recuperación varía entre individuos y sitios corporales."),
                                        simpson           = tags$p("Con Simpson el efecto de los antibióticos se traduce en un aumento de dominancia — pocas especies resistentes colonizan el espacio dejado por las sensibles."),
                                        observed_features = tags$p("La riqueza cae abruptamente con los antibióticos y se recupera gradualmente. Es el índice que muestra la dinámica temporal del impacto de forma más clara."),
                                        faith_pd          = tags$p("Una caída de Faith's PD durante el tratamiento indica pérdida de linajes filogenéticos enteros — un impacto más profundo y más lento de recuperar que la riqueza taxonómica.")))
    
    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretación — ", nom, " · ", ylb, " · ", ltit)),
        tags$div(style="margin-top:0.75rem;", texto)
    )
  })
  
  # ── Taxonómica ────────────────────────────────────────────────────────
  output$p_tax <- renderPlot({
    req(input$sit_tax)
    
    dat <- tax_data %>%
      filter(body_site %in% input$sit_tax) %>%
      mutate(sitio_f = factor(sitio_lbl[as.character(body_site)],
                              levels=sitio_lbl[input$sit_tax]),
             phylum  = factor(phylum, levels=names(pal_phylum)))
    
    if (input$orient_tax == "v") {
      ggplot(dat, aes(x=sitio_f, y=abundancia, fill=phylum)) +
        geom_col(width=0.68) +
        scale_fill_manual(values=pal_phylum, name="Phylum") +
        scale_y_continuous(labels=label_percent(scale=1),
                           expand=expansion(mult=c(0,0.02))) +
        labs(title="Composición taxonómica por sitio corporal",
             subtitle="Abundancia relativa media · Silva 138",
             x=NULL, y="Abundancia relativa (%)") +
        tema_app + theme(legend.position="right")
    } else {
      ggplot(dat, aes(y=sitio_f, x=abundancia, fill=phylum)) +
        geom_col(width=0.68) +
        scale_fill_manual(values=pal_phylum, name="Phylum") +
        scale_x_continuous(labels=label_percent(scale=1),
                           expand=expansion(mult=c(0,0.02))) +
        labs(title="Composición taxonómica por sitio corporal",
             subtitle="Abundancia relativa media · Silva 138",
             y=NULL, x="Abundancia relativa (%)") +
        tema_app + theme(legend.position="bottom")
    }
  })
}

shinyApp(ui, server)