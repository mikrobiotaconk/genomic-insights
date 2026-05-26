# Genomic Insights — Arquitectura Blog + Apps

## Estructura del proyecto

Dos repositorios/carpetas separados que trabajan juntos:

```
genomic-insights/               ← BLOG ESTÁTICO (Quarto website)
├── _quarto.yml                 configuración global
├── index.qmd                   página de inicio + listado automático
├── about.qmd                   página profesional
├── styles.css                  estilos editoriales compartidos
├── assets/
│   ├── avatar.jpg              tu foto (añadir)
│   └── favicon.ico             (añadir)
└── posts/
    ├── _template/              plantilla para nuevos posts
    │   └── index.qmd
    └── moving-pictures/        Post 1
        ├── index.qmd           texto + figura estática + iframe de la app
        └── thumbnail.png       imagen de portada (añadir, ~800×400 px)

genomic-insights-apps/          ← APPS SHINY (Posit Connect)
├── _template/                  plantilla para nuevas apps
│   └── app.R
└── moving-pictures/            App del Post 1
    ├── app.R                   app Shiny completa (4 figuras interactivas)
    └── data/                   ← carpeta para tus datos exportados de QIIME2
        └── (tus .tsv aquí)
```

---

## Flujo de trabajo para un nuevo post

```
1. Analiza en QIIME2 / R
2. Exporta los artefactos a .tsv
         ↓
3. Crea la app Shiny:
   cp -r genomic-insights-apps/_template genomic-insights-apps/nuevo-analisis
   # edita app.R, carga tus datos, diseña las figuras interactivas
         ↓
4. Publica la app en Posit Connect → obtienes una URL
         ↓
5. Crea el post estático:
   cp -r genomic-insights/posts/_template genomic-insights/posts/nuevo-analisis
   # edita index.qmd, pega la URL de la app en APP_URL
         ↓
6. Publica el blog → el nuevo post aparece automáticamente en el índice
```

---

## Setup inicial

### Requisitos R

```r
install.packages(c(
  "shiny", "bslib",                          # UI de las apps
  "ggplot2", "dplyr", "tidyr", "scales",     # visualización
  "vegan",                                   # diversidad microbiana
  "rsconnect",                               # publicar en Posit Connect
  "quarto"                                   # renderizar el blog
))

# Para análisis reales con QIIME2:
# BiocManager::install(c("phyloseq", "DESeq2", "microbiome"))
```

### Versión mínima de Quarto: 1.4

```bash
quarto --version
# Si es < 1.4: descarga en https://quarto.org/docs/get-started/
```

---

## Desarrollo local

```bash
# Blog estático — previsualizar con hot-reload
cd genomic-insights
quarto preview

# App Shiny — lanzar localmente
cd genomic-insights-apps/moving-pictures
Rscript -e "shiny::runApp()"
```

---

## Publicar en Posit Connect

### 1. Configurar cuenta (una sola vez)

```r
library(rsconnect)
rsconnect::setAccountInfo(
  name   = "tu-usuario",   # tu usuario en posit.cloud
  token  = "TU_TOKEN",     # posit.cloud > Account > Tokens > Add Token
  secret = "TU_SECRET"
)
```

### 2. Publicar la Shiny app PRIMERO

```r
rsconnect::deployApp(
  appDir  = "genomic-insights-apps/moving-pictures",
  appName = "moving-pictures-app",
  account = "tu-usuario"
)
# → Anota la URL: https://tu-usuario.posit.cloud/moving-pictures-app
```

### 3. Pegar la URL en el post

Abre `genomic-insights/posts/moving-pictures/index.qmd` y edita:

```r
APP_URL <- "https://tu-usuario.posit.cloud/moving-pictures-app"
```

### 4. Publicar el blog

```r
rsconnect::deployApp(
  appDir  = "genomic-insights",
  appName = "genomic-insights",
  account = "tu-usuario"
)
```

O desde la terminal de RStudio, dentro de `genomic-insights/`:

```bash
quarto publish connect
```

---

## Conectar tus datos reales de QIIME2

### Exportar desde QIIME2

```bash
# Diversidad alfa
qiime tools export \
  --input-path  diversity-metrics/shannon_vector.qza \
  --output-path genomic-insights-apps/moving-pictures/data/

# PCoA Bray-Curtis
qiime tools export \
  --input-path  diversity-metrics/bray_curtis_pcoa_results.qza \
  --output-path genomic-insights-apps/moving-pictures/data/

# Tabla taxonómica colapsada a nivel de phylum
qiime taxa collapse \
  --i-table feature-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table phylum-table.qza
qiime tools export \
  --input-path  phylum-table.qza \
  --output-path genomic-insights-apps/moving-pictures/data/
```

### Cargar en app.R

```r
# En el bloque de DATOS de app.R, sustituye los datos simulados por:

alpha <- readr::read_tsv("data/alpha-diversity.tsv", skip=1,
           col_names=c("sample_id","shannon")) %>%
         left_join(metadata, by="sample_id")   # metadata: tu archivo .tsv de QIIME2

pcoa_bc <- readr::read_delim("data/ordination.txt",
             skip=9, delim="\t", col_names=FALSE) %>%
           select(sample_id=1, PC1=2, PC2=3)
```

---

## Personalización

| Qué | Dónde |
|-----|-------|
| Nombre del blog | `_quarto.yml` → `website: title` |
| URLs de redes | `_quarto.yml` → `navbar: right` y `about.qmd` |
| Colores principales | `styles.css` → `:root { --verde: ... }` |
| Tu foto | `assets/avatar.jpg` |
| Tu bio | `about.qmd` |
| Imagen de portada de cada post | `posts/[post]/thumbnail.png` |

---

## Añadir un nuevo post (resumen)

```bash
# Blog
cp -r genomic-insights/posts/_template genomic-insights/posts/mi-nuevo-post

# App
cp -r genomic-insights-apps/_template genomic-insights-apps/mi-nuevo-post
```

Edita ambos archivos, publica la app, copia la URL al post, publica el blog.
El nuevo post aparece automáticamente en el índice.
