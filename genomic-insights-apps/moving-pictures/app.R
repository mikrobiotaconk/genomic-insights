# ════════════════════════════════════════════════════════════════════════
# Moving Pictures Microbiome Explorer — app_en.R
# English version
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

# ── Palettes ──────────────────────────────────────────────────────────────
pal_body_site <- c("gut"="#1D9E75","left palm"="#EF9F27",
                   "right palm"="#D4537E","tongue"="#378ADD")
pal_subject   <- c("subject-1"="#1D9E75","subject-2"="#378ADD")
pal_antibiotic<- c("Yes"="#D4537E","No"="#1D9E75")

get_pal <- function(var) switch(var,
                                body_site  = pal_body_site,
                                subject    = pal_subject,
                                antibiotic = pal_antibiotic)

get_legend_title <- function(var) switch(var,
                                         body_site  = "Body site",
                                         subject    = "Subject",
                                         antibiotic = "Antibiotic use")

get_labels <- function(var) switch(var,
                                   body_site  = c("gut"="Gut","left palm"="Left palm",
                                                  "right palm"="Right palm","tongue"="Tongue"),
                                   subject    = c("subject-1"="Subject 1","subject-2"="Subject 2"),
                                   antibiotic = c("Yes"="With antibiotics","No"="Without antibiotics"))

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
# DATA
# ════════════════════════════════════════════════════════════════════════
physeq <- readRDS("data/physeq.rds")

# ── Metadata ──────────────────────────────────────────────────────────────
meta <- sample_data(physeq) %>%
  data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  rename(body_site  = body.site,
         antibiotic = reported.antibiotic.usage,
         days       = days.since.experiment.start) %>%
  mutate(body_site = factor(body_site,
                            levels = c("gut","left palm","right palm","tongue")),
         subject   = factor(subject))

# ── Alpha diversity ───────────────────────────────────────────────────────
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

# ── Beta diversity ────────────────────────────────────────────────────────
ord_bc <- ordinate(physeq, method="PCoA", distance="bray")
ord_wu <- ordinate(physeq, method="PCoA", distance="wunifrac")
ord_uu <- ordinate(physeq, method="PCoA", distance="unifrac")

extract_pcoa <- function(ord, metric) {
  coords <- as.data.frame(ord$vectors[,1:2])
  colnames(coords) <- c("Ax1","Ax2")
  pct <- round(ord$values$Relative_eig[1:2]*100,1)
  coords %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(meta, by="sample_id") %>%
    mutate(metric=metric, pct1=pct[1], pct2=pct[2])
}

pcoa_bc <- extract_pcoa(ord_bc, "bc")
pcoa_wu <- extract_pcoa(ord_wu, "wu")
pcoa_uu <- extract_pcoa(ord_uu, "uu")

# ── Taxonomy ──────────────────────────────────────────────────────────────
ps_phylum <- tax_glom(physeq, taxrank="Phylum")

tax_data <- psmelt(ps_phylum) %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance/sum(Abundance)*100) %>%
  ungroup() %>%
  rename(sample_id = Sample) %>%
  mutate(body_site = factor(body.site,
                            levels=c("gut","left palm","right palm","tongue"))) %>%
  group_by(body_site, Phylum) %>%
  summarise(abundance=mean(rel_abund), .groups="drop") %>%
  group_by(body_site) %>%
  mutate(phylum_clean = ifelse(
    Phylum %in% names(sort(tapply(abundance,Phylum,mean),
                           decreasing=TRUE))[1:5], Phylum, "Other")) %>%
  group_by(body_site, phylum_clean) %>%
  summarise(abundance=sum(abundance), .groups="drop") %>%
  rename(phylum=phylum_clean) %>%
  mutate(body_site=factor(body_site,
                          levels=c("gut","left palm","right palm","tongue")))

phyla_unique <- unique(tax_data$phylum)
pal_phylum   <- setNames(
  c("#1D9E75","#378ADD","#FAC775","#F09595","#5DCAA5","#B4B2A9",
    "#9B59B6","#E67E22","#1ABC9C","#E74C3C")[seq_len(length(phyla_unique))],
  phyla_unique)

# ── Helpers ───────────────────────────────────────────────────────────────
site_lbl <- c("gut"="Gut","left palm"="Left palm",
              "right palm"="Right palm","tongue"="Tongue")

label_index <- function(x) switch(x,
                                   shannon="Shannon (H')",
                                   simpson="Simpson (1-D)",
                                   observed_features="Observed richness (ASVs)",
                                   faith_pd="Faith's PD")

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

  # ── Tab 1: Alpha diversity ────────────────────────────────────────────
  nav_panel("Alpha diversity",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controls",
                selectInput("ind_alpha", "Diversity index:",
                            choices  = c("Shannon"="shannon","Simpson (1-D)"="simpson",
                                         "Observed richness"="observed_features","Faith PD"="faith_pd"),
                            selected = "shannon"),
                selectInput("color_alpha", "Colour by:",
                            choices  = c("Body site"="body_site","Subject"="subject",
                                         "Antibiotic use"="antibiotic"),
                            selected = "body_site"),
                checkboxGroupInput("sit_alpha", "Sites:",
                                   choices  = c("Gut"="gut","Tongue"="tongue",
                                                "Left palm"="left palm","Right palm"="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                checkboxGroupInput("suj_alpha", "Subjects:",
                                   choices  = c("Subject 1"="subject-1","Subject 2"="subject-2"),
                                   selected = c("subject-1","subject-2")),
                radioButtons("geom_alpha", "Plot type:",
                             choices  = c("Boxplot"="box","Violin"="violin","Dots + mean"="dots"),
                             selected = "box"),
                checkboxInput("puntos_alpha", "Show individual points", value=TRUE)
              ),
              card(
                full_screen = TRUE,
                card_header("Alpha diversity by body site"),
                plotOutput("p_alpha", height="460px"),
                card_footer(class="text-muted small fst-italic",
                            "Rarefaction: 1,080 reads.")
              ),
              card(
                card_header("📊 Statistical interpretation"),
                uiOutput("texto_alpha")
              )
            )
  ),

  # ── Tab 2: Beta diversity (PCoA) ──────────────────────────────────────
  nav_panel("Beta diversity (PCoA)",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controls",
                selectInput("metrica", "Distance metric:",
                            choices  = c("Bray-Curtis"="bc","Weighted UniFrac"="wu",
                                         "Unweighted UniFrac"="uu"),
                            selected = "bc"),
                selectInput("color_beta", "Colour by:",
                            choices  = c("Body site"="body_site","Subject"="subject",
                                         "Antibiotic use"="antibiotic"),
                            selected = "body_site"),
                checkboxGroupInput("sit_beta", "Sites:",
                                   choices  = c("Gut"="gut","Tongue"="tongue",
                                                "Left palm"="left palm","Right palm"="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                checkboxGroupInput("suj_beta", "Subjects:",
                                   choices  = c("Subject 1"="subject-1","Subject 2"="subject-2"),
                                   selected = c("subject-1","subject-2")),
                checkboxInput("elipses",   "95% confidence ellipses", value=TRUE),
                checkboxInput("etiquetas", "Label samples",           value=FALSE)
              ),
              card(
                full_screen = TRUE,
                card_header("PCoA ordination"),
                plotOutput("p_pcoa", height="460px"),
                card_footer(class="text-muted small fst-italic",
                            "PERMANOVA: R² = 0.68, p = 0.001 (body site).")
              ),
              card(
                card_header("📊 Statistical interpretation"),
                uiOutput("texto_beta")
              )
            )
  ),

  # ── Tab 3: Temporal evolution ─────────────────────────────────────────
  nav_panel("Temporal evolution",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controls",
                selectInput("sit_temp", "Body site:",
                            choices  = c("Gut"="gut","Tongue"="tongue",
                                         "Left palm"="left palm","Right palm"="right palm"),
                            selected = "gut"),
                selectInput("ind_temp", "Index:",
                            choices  = c("Shannon"="shannon","Simpson (1-D)"="simpson",
                                         "Observed richness"="observed_features","Faith PD"="faith_pd"),
                            selected = "shannon"),
                selectInput("color_temp", "Colour by:",
                            choices  = c("Subject"="subject","Antibiotic use"="antibiotic"),
                            selected = "subject"),
                checkboxInput("loess",  "Trend line (loess)", value=TRUE),
                checkboxInput("ribbon", "Confidence interval", value=TRUE),
                hr(),
                helpText("Each subject maintains a personal diversity 'set point'.")
              ),
              card(
                full_screen = TRUE,
                card_header("Temporal dynamics of diversity"),
                plotOutput("p_temp", height="460px"),
                card_footer(class="text-muted small fst-italic",
                            "X axis: days since experiment start.")
              ),
              card(
                card_header("📊 Statistical interpretation"),
                uiOutput("texto_temp")
              )
            )
  ),

  # ── Tab 4: Taxonomic composition ─────────────────────────────────────
  nav_panel("Taxonomic composition",
            layout_sidebar(
              sidebar = sidebar(
                width = 220, title = "Controls",
                checkboxGroupInput("sit_tax", "Sites:",
                                   choices  = c("Gut"="gut","Tongue"="tongue",
                                                "Left palm"="left palm","Right palm"="right palm"),
                                   selected = c("gut","tongue","left palm","right palm")),
                radioButtons("orient_tax", "Orientation:",
                             choices = c("Vertical"="v","Horizontal"="h"), selected="v"),
                hr(),
                helpText("Top 5 phyla per site; remainder grouped as 'Other'.")
              ),
              card(
                full_screen = TRUE,
                card_header("Relative phylum abundance by body site"),
                plotOutput("p_tax", height="460px"),
                card_footer(class="text-muted small fst-italic",
                            "Mean relative abundance (both subjects). Classification: Silva 138.")
              )
            )
  )
)

# ════════════════════════════════════════════════════════════════════════
# SERVER
# ════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Alpha ──────────────────────────────────────────────────────────────
  output$p_alpha <- renderPlot({
    req(input$suj_alpha, input$sit_alpha)

    dat <- alpha_df %>%
      filter(subject %in% input$suj_alpha, body_site %in% input$sit_alpha) %>%
      mutate(sitio_f   = factor(site_lbl[as.character(body_site)],
                                levels=site_lbl[input$sit_alpha]),
             color_var = as.character(.data[[input$color_alpha]]))

    y <- sym(input$ind_alpha); ylb <- label_index(input$ind_alpha)
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
      labs(title=paste("Alpha diversity —", ylb),
           subtitle=paste("Coloured by:", ltit), x=NULL, y=ylb) +
      tema_app + theme(legend.position="top")
  })

  # ── Alpha interpretation text ─────────────────────────────────────────
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
        Comparison  = dunn$comparisons,
        Z           = round(dunn$Z, 3),
        p_adjusted  = format.pval(dunn$P.adjusted, digits=3, eps=0.001),
        significant = ifelse(dunn$P.adjusted < 0.05, "✓", "—")
      )
      rows <- apply(dunn_df, 1, function(r) {
        color_sig <- if (r["significant"]=="✓") "#085041" else "#aaa"
        tags$tr(
          tags$td(r["Comparison"],  style="padding:4px 8px; font-size:0.82rem;"),
          tags$td(r["Z"],           style="padding:4px 8px; font-size:0.82rem; text-align:center;"),
          tags$td(r["p_adjusted"],  style="padding:4px 8px; font-size:0.82rem; text-align:center;"),
          tags$td(r["significant"], style=paste0("padding:4px 8px; font-size:0.82rem;
                  text-align:center; font-weight:600; color:",color_sig))
        )
      })
      dunn_html <- tagList(
        tags$p(style="font-size:0.78rem; font-weight:500; color:#085041; margin:0.75rem 0 0.25rem;",
               "Post-hoc comparisons (Dunn · Bonferroni correction):"),
        tags$table(style="width:100%; border-collapse:collapse;",
                   tags$thead(tags$tr(style="border-bottom:1px solid #9FE1CB;",
                                      tags$th("Comparison",  style="padding:4px 8px; font-size:0.78rem; text-align:left;"),
                                      tags$th("Z",           style="padding:4px 8px; font-size:0.78rem; text-align:center;"),
                                      tags$th("Adj. p",      style="padding:4px 8px; font-size:0.78rem; text-align:center;"),
                                      tags$th("Sig.",        style="padding:4px 8px; font-size:0.78rem; text-align:center;")
                   )),
                   tags$tbody(rows)
        )
      )
    }

    texto <- switch(input$color_alpha,
                    body_site = switch(input$ind_alpha,
                                       shannon           = tags$p("Shannon weights both richness and evenness. The oral cavity leads (H'~3.5), followed by the gut (~3.0) and skin (~2.0)."),
                                       simpson           = tags$p("Simpson gives more weight to dominant species. Patterns hold but differences are less pronounced than with Shannon."),
                                       observed_features = tags$p("Raw ASV richness follows the same gradient but is more sensitive to sequencing depth."),
                                       faith_pd          = tags$p("Faith's PD incorporates evolutionary distance between lineages — it reveals phylogenetic diversity, not just taxonomic."),
                                       tags$p("Select an index.")),
                    subject = switch(input$ind_alpha,
                                     shannon           = tags$p("Variation between subjects is smaller than between sites. Each individual maintains a stable personal diversity 'set point' over time."),
                                     simpson           = tags$p("With Simpson, inter-subject differences are even less pronounced — dominance of major species is similar across individuals."),
                                     observed_features = tags$p("Raw richness varies more between subjects than weighted indices, being more sensitive to differences in specific composition."),
                                     faith_pd          = tags$p("Phylogenetic diversity between subjects is surprisingly similar despite differing in taxonomic composition."),
                                     tags$p("Select an index.")),
                    antibiotic = switch(input$ind_alpha,
                                        shannon           = tags$p("Antibiotics clearly reduce Shannon, especially in the gut. Recovery is visible in subsequent samples."),
                                        simpson           = tags$p("With Simpson the effect is more visible: dominance of a few resistant species increases during treatment."),
                                        observed_features = tags$p("Richness drops abruptly — entire ASVs are lost. This is the most sensitive index for detecting the initial impact."),
                                        faith_pd          = tags$p("A drop in Faith's PD indicates loss of entire evolutionary lineages — harder to recover than taxonomic richness."),
                                        tags$p("Select an index.")),
                    tags$p("Select a variable."))

    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretation — ", get_legend_title(input$color_alpha),
                           " · ", label_index(input$ind_alpha))),
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
      mutate(sitio_f   = factor(site_lbl[as.character(body_site)],
                                levels=site_lbl[input$sit_beta]),
             color_var = as.character(.data[[input$color_beta]]))

    pct1 <- dat$pct1[1]; pct2 <- dat$pct2[1]
    nom  <- switch(input$metrica, bc="Bray-Curtis",
                   wu="Weighted UniFrac", uu="Unweighted UniFrac")
    pal  <- get_pal(input$color_beta)
    lbl  <- get_labels(input$color_beta)
    ltit <- get_legend_title(input$color_beta)

    p <- ggplot(dat, aes(x=Ax1, y=Ax2, color=color_var, fill=color_var)) +
      geom_point(size=4, alpha=0.85) +
      scale_color_manual(values=pal, name=ltit, labels=lbl) +
      scale_fill_manual( values=pal, name=ltit, labels=lbl, guide="none") +
      labs(title=paste0("PCoA (", nom, ")"),
           subtitle=paste("Coloured by:", ltit),
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

  # ── Beta interpretation text ──────────────────────────────────────────
  output$texto_beta <- renderUI({
    req(input$metrica, input$color_beta)

    nom  <- switch(input$metrica, bc="Bray-Curtis",
                   wu="Weighted UniFrac", uu="Unweighted UniFrac")
    ltit <- get_legend_title(input$color_beta)

    texto <- switch(input$color_beta,
                    body_site = switch(input$metrica,
                                       bc = tags$p("Bray-Curtis clearly separates the four body sites. Site explains 68% of total variance (PERMANOVA R²=0.68, p=0.001) — far more than the individual subject (R²=0.09)."),
                                       wu = tags$p("Weighted UniFrac incorporates relative abundance and phylogenetic distance. Separation between sites holds but reflects deeper evolutionary differences between communities."),
                                       uu = tags$p("Unweighted UniFrac is more sensitive to rare species. It can reveal inter-site differences that Bray-Curtis misses, as it is not dominated by the most abundant species.")),
                    subject = switch(input$metrica,
                                     bc = tags$p("With Bray-Curtis, samples from the same subject tend to cluster within each body site, but subject separation is much weaker than site separation (R²=0.09 vs 0.68)."),
                                     wu = tags$p("Weighted UniFrac can reveal more subtle inter-subject differences at the phylogenetic level, although the dominant pattern remains body site."),
                                     uu = tags$p("Unweighted UniFrac sometimes amplifies inter-subject differences by giving more weight to rare species that may be individual-specific.")),
                    antibiotic = switch(input$metrica,
                                        bc = tags$p("Antibiotic use can shift samples in PCoA space, especially in the gut. Samples collected during treatment tend to move away from their usual site cluster."),
                                        wu = tags$p("Weighted UniFrac captures the phylogenetic impact of antibiotics well — loss of entire lineages shifts samples further than simple abundance reduction."),
                                        uu = tags$p("Unweighted UniFrac is especially sensitive to the antibiotic effect on rare species, which may disappear completely during treatment.")))

    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretation — ", ltit, " · ", nom)),
        tags$div(style="margin-top:0.75rem;", texto)
    )
  })

  # ── Temporal ──────────────────────────────────────────────────────────
  output$p_temp <- renderPlot({
    dat <- alpha_df %>%
      filter(body_site == input$sit_temp) %>%
      mutate(color_var = as.character(.data[[input$color_temp]]))

    y <- sym(input$ind_temp); ylb <- label_index(input$ind_temp)
    nom  <- site_lbl[input$sit_temp]
    pal  <- get_pal(input$color_temp)
    lbl  <- get_labels(input$color_temp)
    ltit <- get_legend_title(input$color_temp)

    p <- ggplot(dat, aes(x=days, y=!!y, color=color_var, group=color_var)) +
      geom_line(linewidth=0.9, alpha=0.75) +
      geom_point(size=2.5, alpha=0.85) +
      scale_color_manual(values=pal, name=ltit, labels=lbl) +
      labs(title=paste0("Temporal variation — ", nom),
           subtitle=paste0(ylb, " · coloured by: ", ltit),
           x="Days since experiment start", y=ylb) +
      tema_app + theme(legend.position="top")

    if (input$loess)
      p <- p + geom_smooth(aes(group=color_var), method="loess",
                           se=input$ribbon, linewidth=0.55,
                           alpha=0.12, linetype="dashed")
    p
  })

  # ── Temporal interpretation text ──────────────────────────────────────
  output$texto_temp <- renderUI({
    req(input$sit_temp, input$ind_temp, input$color_temp)

    nom  <- site_lbl[input$sit_temp]
    ylb  <- label_index(input$ind_temp)
    ltit <- get_legend_title(input$color_temp)

    texto <- switch(input$color_temp,
                    subject = switch(input$ind_temp,
                                     shannon           = tags$p("Each subject maintains a relatively stable personal Shannon 'set point'. Fluctuations are moderate and the system tends to recover its baseline after perturbations."),
                                     simpson           = tags$p("With Simpson, individual stability is even more pronounced — dominance of the main species is very consistent over time within each individual."),
                                     observed_features = tags$p("Observed richness fluctuates more than weighted indices, being more sensitive to transient events such as dietary changes or mild infections."),
                                     faith_pd          = tags$p("Phylogenetic diversity shows the greatest temporal stability — the evolutionary lineages present change little even as individual species shift.")),
                    antibiotic = switch(input$ind_temp,
                                        shannon           = tags$p("Antibiotic treatment causes visible short-term drops in Shannon along the time axis. Recovery speed varies between individuals and body sites."),
                                        simpson           = tags$p("With Simpson, the antibiotic effect translates into increased dominance — a few resistant species colonise the space left by sensitive ones."),
                                        observed_features = tags$p("Richness drops abruptly with antibiotics and recovers gradually. This is the index that shows the temporal impact dynamics most clearly."),
                                        faith_pd          = tags$p("A drop in Faith's PD during treatment indicates loss of entire phylogenetic lineages — a deeper and slower-to-recover impact than taxonomic richness.")))

    div(style="background:#E1F5EE; border-left:4px solid #1D9E75;
               padding:1rem 1.25rem; border-radius:0 8px 8px 0;
               margin:0.5rem; font-size:0.92rem; line-height:1.75;",
        tags$strong(style="font-size:0.75rem; text-transform:uppercase;
                          letter-spacing:0.08em; color:#085041;",
                    paste0("Interpretation — ", nom, " · ", ylb, " · ", ltit)),
        tags$div(style="margin-top:0.75rem;", texto)
    )
  })

  # ── Taxonomic composition ─────────────────────────────────────────────
  output$p_tax <- renderPlot({
    req(input$sit_tax)

    dat <- tax_data %>%
      filter(body_site %in% input$sit_tax) %>%
      mutate(sitio_f = factor(site_lbl[as.character(body_site)],
                              levels=site_lbl[input$sit_tax]),
             phylum  = factor(phylum, levels=names(pal_phylum)))

    if (input$orient_tax == "v") {
      ggplot(dat, aes(x=sitio_f, y=abundance, fill=phylum)) +
        geom_col(width=0.68) +
        scale_fill_manual(values=pal_phylum, name="Phylum") +
        scale_y_continuous(labels=label_percent(scale=1),
                           expand=expansion(mult=c(0,0.02))) +
        labs(title="Taxonomic composition by body site",
             subtitle="Mean relative abundance · Silva 138",
             x=NULL, y="Relative abundance (%)") +
        tema_app + theme(legend.position="right")
    } else {
      ggplot(dat, aes(y=sitio_f, x=abundance, fill=phylum)) +
        geom_col(width=0.68) +
        scale_fill_manual(values=pal_phylum, name="Phylum") +
        scale_x_continuous(labels=label_percent(scale=1),
                           expand=expansion(mult=c(0,0.02))) +
        labs(title="Taxonomic composition by body site",
             subtitle="Mean relative abundance · Silva 138",
             y=NULL, x="Relative abundance (%)") +
        tema_app + theme(legend.position="bottom")
    }
  })
}

shinyApp(ui, server)
