# =============================================================================
# Cartographie des Myriapodes et Isopodes de France
# -----------------------------------------------------------------------------
# Sources de données : GBIF (API), CSV (export PostgreSQL)
# Référentiel taxonomique principal : GBIF Backbone Taxonomy
# Photos : Wikimedia Commons (CC), avec repli sur les médias GBIF
#
# IMPORTANT :
# - GBIF est utilisé pour la recherche taxonomique et les occurrences GBIF.
# - Le CSV est rapproché des taxons par nom scientifique (nom GBIF ou
#   nom canonique). Les identifiants cd_ref du CSV restent présents mais ne
#   sont pas utilisés pour sélectionner le taxon.
# - OpenObs n'est pas interrogé dans cette version : son paramètre taxon
#   attend un identifiant TAXREF, qui ne peut pas être remplacé par usageKey.
#   Pour réactiver OpenObs, il faudra une table de correspondance GBIF/TAXREF.
#
# Lancement : shiny::runApp("app.R")
# Dépendances :
# install.packages(c("shiny", "leaflet", "httr", "jsonlite", "dplyr"))
#
# CSV attendu : data/synthese_observations.csv (séparateur ";")
# Colonnes minimales : nom_valide ; cd_ref ; x_centroid_4326 ; y_centroid_4326
# Colonnes facultatives : observateurs ; date_debut ; nombre_max ; sexe ;
#                         nom_vernaculaire
# =============================================================================

library(shiny)
library(leaflet)
library(httr)
library(jsonlite)
library(dplyr)

# Valeur de repli pour NULL, vecteur vide ou NA
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
}

# ----------------------------- Constantes ------------------------------------

CSV_PATH <- "C:/Users/Thomas/OneDrive - ANA - Conservatoire d'espaces naturels Ariege/Invertébrés/Myriapodes/MyriaFrance/shiny/synthese_observations_2026-11-13T10_58_39.118Z.csv"
MAX_OCC_GBIF <- 3000
GBIF_PAGE_SIZE <- 300
MAX_ZOOM <- 19

PALETTE <- c(
  "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
  "#00a2c7", "#f032e6", "#bcf60c", "#fabebe", "#008080",
  "#9a6324", "#800000", "#808000", "#000075", "#e6beff"
)

# Formes : carré = GBIF, cercle = OpenObs (réservé), triangle = CSV
SHAPES <- c(gbif = "square", inpn = "circle", csv = "triangle")

# ----------------------------- Utilitaires -----------------------------------

get_json <- function(url, query = list()) {
  res <- tryCatch(
    httr::GET(
      url,
      query = query,
      httr::timeout(30),
      httr::user_agent("shiny-myriapodes-app/2.0")
    ),
    error = function(e) {
      message("Erreur réseau pour ", url, " : ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res)) return(NULL)

  if (httr::status_code(res) != 200) {
    message(
      "Réponse HTTP ", httr::status_code(res),
      " pour ", url
    )
    return(NULL)
  }

  tryCatch(
    jsonlite::fromJSON(
      httr::content(res, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    ),
    error = function(e) {
      message("Erreur de lecture JSON : ", conditionMessage(e))
      NULL
    }
  )
}

safe_text <- function(x, default = "") {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) default else as.character(x[1])
}

# ----------------------------- GBIF TAXONOMY ---------------------------------

# Autocomplétion via le référentiel GBIF Backbone.
# L'API /species/suggest renvoie une liste de taxons correspondant à la saisie.
# Vérifie qu'un taxon appartient aux Chilopoda, Diplopoda ou à Isopoda
# en examinant sa lignée dans le GBIF Backbone.
gbif_is_target_group <- function(usageKey) {
  if (is.null(usageKey) || length(usageKey) == 0 || is.na(usageKey)) {
    return(FALSE)
  }

  self <- get_json(sprintf("https://api.gbif.org/v1/species/%s", usageKey))
  js <- get_json(
    sprintf("https://api.gbif.org/v1/species/%s/parents", usageKey)
  )
  if (is.null(js) && is.null(self)) return(FALSE)

  ancestors <- if (is.list(js) && !is.null(js$results)) js$results else js
  if (is.null(ancestors)) ancestors <- list()
  if (!is.list(ancestors)) ancestors <- list(ancestors)
  if (!is.null(self)) ancestors <- c(list(self), ancestors)

  ancestor_names <- vapply(ancestors, function(x) {
    tolower(trimws(safe_text(x$scientificName %||% x$canonicalName)))
  }, character(1))

  any(ancestor_names %in% c("chilopoda", "diplopoda", "isopoda"))
}

# Récupère récursivement les noms des taxons descendants d'un genre.
# Ces noms servent à rapprocher les observations du CSV.
gbif_descendant_names <- function(usageKey, max_nodes = 5000) {
  found <- character()
  visited <- character()

  walk <- function(parent_key, depth = 0L) {
    if (length(visited) >= max_nodes || depth > 8L) return(invisible(NULL))
    parent_key <- as.character(parent_key)
    if (parent_key %in% visited) return(invisible(NULL))
    visited <<- c(visited, parent_key)

    js <- get_json(
      sprintf("https://api.gbif.org/v1/species/%s/children", parent_key),
      list(limit = 1000, offset = 0)
    )
    if (is.null(js)) return(invisible(NULL))

    children <- if (is.list(js) && !is.null(js$results)) js$results else js
    if (!is.list(children) || length(children) == 0) return(invisible(NULL))

    # Une réponse GBIF peut être une liste d'objets ou un objet unique.
    if (!is.null(children$key)) children <- list(children)

    for (child in children) {
      if (length(visited) >= max_nodes) break
      child_key <- child$key %||% NULL
      child_name <- safe_text(child$scientificName)
      child_canonical <- safe_text(child$canonicalName)

      found <<- c(found, child_name, child_canonical)

      # Explore aussi les descendants d'une espèce afin d'inclure
      # d'éventuels rangs infra-spécifiques.
      if (!is.null(child_key)) {
        walk(child_key, depth + 1L)
      }
    }
    invisible(NULL)
  }

  walk(usageKey)
  unique(tolower(trimws(found[nzchar(found)])))
}

gbif_suggest <- function(q) {
  q <- trimws(as.character(q %||% ""))

  if (length(q) == 0 || is.na(q) || nchar(q) < 3) return(NULL)

  js <- get_json(
    "https://api.gbif.org/v1/species/suggest",
    list(q = q, limit = 30)
  )

  if (is.null(js) || length(js) == 0) {
    match_js <- get_json(
      "https://api.gbif.org/v1/species/match",
      list(name = q, strict = "false")
    )
    if (is.null(match_js) || is.null(match_js$key)) return(NULL)
    js <- list(match_js)
  }

  rows <- lapply(js, function(t) {
    key <- t$key %||% NA
    accepted_key <- t$acceptedUsageKey %||% key

    data.frame(
      usageKey = as.character(key),
      acceptedUsageKey = as.character(accepted_key),
      name = safe_text(t$scientificName),
      canonicalName = safe_text(t$canonicalName),
      rank = safe_text(t$rank),
      status = safe_text(t$status),
      vern = safe_text(t$vernacularName),
      stringsAsFactors = FALSE
    )
  })

  taxa <- dplyr::bind_rows(rows)
  taxa <- taxa[!is.na(taxa$usageKey) & nzchar(taxa$usageKey), , drop = FALSE]
  if (nrow(taxa) == 0) return(NULL)
  taxa <- taxa[!duplicated(taxa$acceptedUsageKey), , drop = FALSE]

  # Filtrage strict : uniquement les taxons dont la lignée contient
  # Chilopoda, Diplopoda ou Isopoda.
  keep <- vapply(taxa$acceptedUsageKey, gbif_is_target_group, logical(1))
  taxa <- taxa[keep, , drop = FALSE]

  if (nrow(taxa) == 0) return(NULL)
  taxa
}

# ----------------------------- GBIF OCCURRENCES -------------------------------

gbif_occurrences <- function(usageKey) {
  if (is.null(usageKey) || length(usageKey) == 0 || is.na(usageKey)) {
    return(NULL)
  }

  all_pts <- list()
  offset <- 0L

  repeat {
    js <- get_json(
      "https://api.gbif.org/v1/occurrence/search",
      list(
        taxonKey = as.character(usageKey),
        hasCoordinate = "true",
        country = "FR",
        limit = GBIF_PAGE_SIZE,
        offset = offset
      )
    )

    if (is.null(js) || is.null(js$results) || length(js$results) == 0) break

    pts <- lapply(js$results, function(o) {
      lat <- suppressWarnings(as.numeric(o$decimalLatitude %||% NA))
      lon <- suppressWarnings(as.numeric(o$decimalLongitude %||% NA))

      if (length(lat) == 0 || length(lon) == 0 ||
          is.na(lat[1]) || is.na(lon[1])) {
        return(NULL)
      }

      data.frame(
        taxon = safe_text(o$scientificName),
        lat = lat[1],
        lon = lon[1],
        source = "gbif",
        cd_ref = "",
        obs = safe_text(o$recordedBy),
        date = substr(safe_text(o$eventDate), 1, 10),
        effectif = safe_text(o$individualCount),
        sexe = "",
        vern = "",
        stringsAsFactors = FALSE
      )
    })

    page_df <- dplyr::bind_rows(pts)
    if (nrow(page_df) > 0) {
      all_pts[[length(all_pts) + 1]] <- page_df
    }

    offset <- offset + GBIF_PAGE_SIZE

    total <- suppressWarnings(as.numeric(js$count %||% 0))
    if (is.na(total) || offset >= total || offset >= MAX_OCC_GBIF) break
  }

  if (length(all_pts) == 0) return(NULL)
  dplyr::bind_rows(all_pts)
}

# ----------------------------- GBIF PHOTO -------------------------------------

gbif_photo <- function(usageKey, max_photos = 12) {
  if (is.null(usageKey) || length(usageKey) == 0 || is.na(usageKey)) {
    return(NULL)
  }

  js <- get_json(
    "https://api.gbif.org/v1/occurrence/search",
    list(
      taxonKey = as.character(usageKey),
      mediaType = "StillImage",
      limit = 100
    )
  )

  if (is.null(js) || is.null(js$results) || length(js$results) == 0) {
    return(NULL)
  }

  photos <- list()
  for (occ in js$results) {
    med <- occ$media %||% NULL
    if (is.null(med) || length(med) == 0) next

    for (m in med) {
      url <- m$identifier %||% NULL
      if (is.null(url) || !nzchar(as.character(url))) next
      photos[[length(photos) + 1]] <- list(
        url = as.character(url),
        credit = safe_text(m$license, "Licence non précisée"),
        source = "GBIF"
      )
      if (length(photos) >= max_photos) break
    }
    if (length(photos) >= max_photos) break
  }

  if (length(photos) == 0) return(NULL)
  urls <- vapply(photos, function(x) x$url, character(1))
  photos[!duplicated(urls)]
}

# ----------------------------- CSV PostgreSQL ---------------------------------

csv_occurrences <- function() {
  if (!file.exists(CSV_PATH)) return(NULL)

  df <- tryCatch(
    read.csv2(
      CSV_PATH,
      fileEncoding = "UTF-8",
      stringsAsFactors = FALSE,
      na.strings = c("", "NA")
    ),
    error = function(e) {
      message("Erreur de lecture du CSV : ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(df)) return(NULL)

  names(df) <- tolower(trimws(names(df)))

  need <- c(
    "nom_valide", "cd_ref",
    "x_centroid_4326", "y_centroid_4326"
  )

  if (!all(need %in% names(df))) {
    warning(
      "CSV : colonnes attendues absentes : ",
      paste(setdiff(need, names(df)), collapse = ", ")
    )
    return(NULL)
  }

  effectif <- if ("nombre_max" %in% names(df)) {
    as.character(df$nombre_max)
  } else {
    rep("", nrow(df))
  }

  data.frame(
    taxon = as.character(df$nom_valide),
    cd_ref = as.character(df$cd_ref),
    lat = suppressWarnings(as.numeric(df$y_centroid_4326)),
    lon = suppressWarnings(as.numeric(df$x_centroid_4326)),
    source = "csv",
    obs = if ("observateurs" %in% names(df)) {
      as.character(df$observateurs)
    } else rep("", nrow(df)),
    date = if ("date_debut" %in% names(df)) {
      as.character(df$date_debut)
    } else rep("", nrow(df)),
    effectif = effectif,
    sexe = if ("sexe" %in% names(df)) {
      as.character(df$sexe)
    } else rep("", nrow(df)),
    vern = if ("nom_vernaculaire" %in% names(df)) {
      as.character(df$nom_vernaculaire)
    } else rep("", nrow(df)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(lat), !is.na(lon),
      lat >= -90, lat <= 90,
      lon >= -180, lon <= 180
    )
}

# ----------------------------- Wikimedia Commons -----------------------------

wikimedia_photo <- function(name, max_photos = 12) {
  if (is.null(name) || !nzchar(name)) return(NULL)

  js <- get_json(
    "https://commons.wikimedia.org/w/api.php",
    list(
      action = "query",
      format = "json",
      origin = "*",
      generator = "search",
      gsrsearch = paste0('"', name, '"'),
      gsrnamespace = 6,
      gsrlimit = max_photos,
      prop = "imageinfo",
      iiprop = "url|extmetadata",
      iiurlwidth = 700
    )
  )

  pages <- js$query$pages %||% NULL
  if (is.null(pages) || length(pages) == 0) return(NULL)

  photos <- list()
  for (pg in pages) {
    ii <- pg$imageinfo[[1]] %||% NULL
    if (is.null(ii)) next

    mime <- safe_text(ii$mime)
    if (!grepl("^image/", mime)) next

    lic <- safe_text(ii$extmetadata$LicenseShortName$value, "Licence non précisée")
    image_url <- ii$thumburl %||% ii$url %||% NULL
    if (is.null(image_url) || !nzchar(as.character(image_url))) next

    photos[[length(photos) + 1]] <- list(
      url = as.character(image_url),
      credit = lic,
      source = "Wikimedia Commons"
    )
    if (length(photos) >= max_photos) break
  }

  if (length(photos) == 0) return(NULL)
  urls <- vapply(photos, function(x) x$url, character(1))
  photos[!duplicated(urls)]
}

# ----------------------------- Icônes SVG ------------------------------------

make_icon <- function(shape, color, size = 22) {
  s <- size
  h <- s / 2

  svg <- switch(
    shape,
    square = sprintf(
      '<rect x="2" y="2" width="%d" height="%d" fill="%s" stroke="white" stroke-width="1.5"/>',
      s - 4, s - 4, color
    ),
    circle = sprintf(
      '<circle cx="%d" cy="%d" r="%d" fill="%s" stroke="white" stroke-width="1.5"/>',
      h, h, h - 2, color
    ),
    triangle = sprintf(
      '<path d="M %d %d L %d %d L %d %d Z" fill="%s" stroke="white" stroke-width="1.5"/>',
      h, 2, s - 2, s - 3, 2, s - 3, color
    )
  )

  svg_full <- sprintf(
    '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d">%s</svg>',
    s, s, svg
  )

  # Le SVG est encodé en base64 : les caractères #, espaces et guillemets
  # ne peuvent ainsi plus casser l'URL de l'icône Leaflet.
  svg_b64 <- jsonlite::base64_enc(charToRaw(enc2utf8(svg_full)))
  paste0("data:image/svg+xml;base64,", svg_b64)
}

# ----------------------------- UI ---------------------------------------------

ui <- fillPage(
  tags$head(
    tags$style(HTML("
      #panel {
        background: rgba(250, 250, 250, 0.97);
        border-radius: 6px;
        box-shadow: 0 1px 6px rgba(0,0,0,0.3);
        padding: 12px;
        box-sizing: border-box;
        max-height: calc(100% - 20px);
        overflow-y: auto;
      }
      #sug_box { position: relative; margin-bottom: 10px; }
      #search { width: 100%; box-sizing: border-box; padding: 8px; }
      .sug-list {
        position: absolute; z-index: 5000; width: 100%;
        background: white; border: 1px solid #ccc; max-height: 260px;
        overflow-y: auto; box-shadow: 0 2px 6px rgba(0,0,0,0.2);
        top: 100%; left: 0;
      }
      .sug-status {
        font-size: 12px; color: #666; padding: 5px 2px;
      }
      .sug-item {
        padding: 6px 8px; cursor: pointer;
        border-bottom: 1px solid #eee; font-size: 13px;
      }
      .sug-item:hover { background: #eef2ff; }
      .taxon-card {
        background: white; border: 1px solid #ddd;
        border-radius: 6px; padding: 10px; margin-bottom: 10px;
      }
      .photo-strip {
        display: flex; gap: 8px; overflow-x: auto; overflow-y: hidden;
        padding: 4px 2px 10px; scroll-snap-type: x mandatory;
        -webkit-overflow-scrolling: touch;
      }
      .photo-slide {
        flex: 0 0 88%; min-width: 0; scroll-snap-align: start;
      }
      .photo-slide img {
        display: block; width: 100%; height: 190px;
        object-fit: contain; background: #f3f3f3; border-radius: 4px;
      }
      .photo-help { text-align: center; margin-top: 2px; }
      .taxon-title {
        font-weight: bold; color: #222; margin-bottom: 4px;
        min-height: 22px; display: flex; align-items: center;
      }
      .taxon-card.collapsed { padding-top: 7px; padding-bottom: 7px; }
      .taxon-card.collapsed .taxon-title { margin-bottom: 0; }
      .taxon-links a { margin-right: 8px; font-size: 12px; }
      .credit { font-size: 11px; color: #888; }
    "))
  ),

  leafletOutput("map", width = "100%", height = "100%"),

  absolutePanel(
    id = "panel",
    top = 10,
    right = 10,
    width = 380,
    draggable = FALSE,

    actionButton(
      "reset_btn",
      "Réinitialiser la carte",
      class = "btn btn-warning btn-sm",
      style = "width:100%; margin-bottom:10px;"
    ),

    div(
      id = "sug_box",
      textInput(
        "search",
        label = NULL,
        placeholder = "Rechercher un taxon (référentiel GBIF)..."
      ),
      uiOutput("sug_list"),
      uiOutput("sug_status")
    ),

    tags$small(
      HTML("<b>Source des observations :</b><br>&#9633; GBIF &nbsp;&nbsp; &#9650; CSV PostgreSQL"),
      style = "display:block; margin-bottom:10px; color:#555;"
    ),

    actionButton(
      "toggle_cards",
      "Réduire les vignettes",
      class = "btn btn-default btn-sm",
      style = "width:100%; margin: 2px 0 8px 0;"
    ),

    uiOutput("taxon_cards")
  )
)

# ----------------------------- Server -----------------------------------------

server <- function(input, output, session) {

  taxa <- reactiveValues(list = list())
  sugs <- reactiveValues(df = NULL)
  sug_message <- reactiveVal("")
  cards_collapsed <- reactiveVal(FALSE)

  observeEvent(input$toggle_cards, {
    cards_collapsed(!cards_collapsed())
    updateActionButton(
      session,
      "toggle_cards",
      label = if (cards_collapsed()) "Afficher les détails des vignettes" else "Réduire les vignettes"
    )
  }, ignoreInit = TRUE)

  # --- Autocomplétion GBIF ---
  observeEvent(
    input$search,
    {
      q <- trimws(input$search %||% "")

      if (nchar(q) < 3) {
        sugs$df <- NULL
        sug_message("")
      } else {
        sug_message("Recherche dans le référentiel GBIF...")
        resultat <- tryCatch(
          gbif_suggest(q),
          error = function(e) {
            message("Erreur autocomplétion GBIF : ", conditionMessage(e))
            NULL
          }
        )
        sugs$df <- resultat

        if (is.null(resultat) || nrow(resultat) == 0) {
          sug_message(
            paste0(
              "Aucune suggestion reçue de GBIF pour « ", q,
              " ». Vérifiez la connexion Internet ou consultez la console R."
            )
          )
        } else {
          sug_message("")
        }
      }
    },
    ignoreInit = TRUE
  )

  output$sug_status <- renderUI({
    msg <- sug_message()
    if (!nzchar(msg)) return(NULL)
    tags$div(class = "sug-status", msg)
  })

  # --- Affichage des suggestions ---
  output$sug_list <- renderUI({
    df <- sugs$df
    q <- trimws(input$search %||% "")

    if (is.null(df) || nrow(df) == 0 || nchar(q) < 3) {
      return(NULL)
    }

    tags$div(
      class = "sug-list",
      lapply(seq_len(nrow(df)), function(i) {
        vern_text <- if (nzchar(df$vern[i])) {
          paste0(" — ", df$vern[i])
        } else {
          ""
        }

        label <- HTML(sprintf(
          "%s <small style='color:#888'>(%s ; %s)%s</small>",
          htmltools::htmlEscape(df$name[i]),
          htmltools::htmlEscape(df$rank[i]),
          htmltools::htmlEscape(df$status[i]),
          htmltools::htmlEscape(vern_text)
        ))

        tags$div(
          class = "sug-item",
          onclick = sprintf(
            "Shiny.setInputValue('sug_pick', %d, {priority:'event'});",
            i
          ),
          label
        )
      })
    )
  })

  # --- Sélection d'un taxon ---
  # IMPORTANT : name, usageKey et les traitements qui en dépendent sont
  # définis à l'intérieur de cet observeEvent. Cela évite l'erreur
  # « objet 'name' introuvable » au démarrage de l'application.
  observeEvent(input$sug_pick, {

    df <- sugs$df
    i <- input$sug_pick

    if (is.null(df) || nrow(df) == 0 ||
        is.null(i) || length(i) == 0 ||
        is.na(i) || i < 1 || i > nrow(df)) {
      return()
    }

    t <- df[i, , drop = FALSE]

    # Informations du taxon sélectionné
    name <- as.character(t$name[1])
    canonicalName <- as.character(t$canonicalName[1])
    usageKey <- suppressWarnings(as.numeric(t$acceptedUsageKey[1]))

    if (is.na(usageKey)) {
      usageKey <- suppressWarnings(as.numeric(t$usageKey[1]))
    }

    if (is.na(usageKey) || !nzchar(name)) {
      showNotification(
        "La clé GBIF de ce taxon est absente ou invalide.",
        type = "error"
      )
      return()
    }

    updateTextInput(session, "search", value = "")
    sugs$df <- NULL
    sug_message("")

    if (name %in% names(taxa$list)) return()

    color <- PALETTE[
      (length(taxa$list) %% length(PALETTE)) + 1
    ]

    # Occurrences GBIF
    occ_gbif <- gbif_occurrences(usageKey)

    # OpenObs est volontairement désactivé :
    # son API attend un identifiant cd_ref TAXREF, pas une clé GBIF.
    occ_inpn <- NULL

    # Occurrences du CSV : pour un genre, inclut tous les noms GBIF
    # descendants (espèces et rangs inférieurs). Pour les autres rangs,
    # conserve le rapprochement exact du nom sélectionné.
    occ_csv <- {
      d <- csv_occurrences()

      if (!is.null(d) && nrow(d) > 0) {
        noms_recherche <- unique(tolower(trimws(c(name, canonicalName))))

        if (identical(toupper(as.character(t$rank[1])), "GENUS")) {
          noms_descendants <- gbif_descendant_names(usageKey)
          noms_recherche <- unique(c(noms_recherche, noms_descendants))
        }

        noms_recherche <- noms_recherche[
          !is.na(noms_recherche) & nzchar(noms_recherche)
        ]
        noms_csv <- tolower(trimws(as.character(d$taxon)))
        sel <- !is.na(noms_csv) & noms_csv %in% noms_recherche
        d[sel, , drop = FALSE]
      } else {
        NULL
      }
    }

    # Plusieurs photographies : on rassemble Wikimedia Commons et GBIF.
    # Les doublons d'URL sont retirés et la galerie est limitée à 12 images.
    photos_wiki <- wikimedia_photo(
      if (nzchar(canonicalName)) canonicalName else name
    )
    photos_gbif <- gbif_photo(usageKey)
    photos <- c(photos_wiki %||% list(), photos_gbif %||% list())
    if (length(photos) > 0) {
      photo_urls <- vapply(photos, function(x) x$url, character(1))
      photos <- photos[!duplicated(photo_urls)]
      photos <- photos[seq_len(min(length(photos), 12))]
    } else {
      photos <- NULL
    }

    # Enregistrement du taxon sélectionné
    taxa$list[[name]] <- list(
      name = name,
      canonicalName = canonicalName,
      rank = as.character(t$rank[1]),
      usageKey = usageKey,
      color = color,
      photos = photos,
      occ = dplyr::bind_rows(
        occ_gbif,
        occ_inpn,
        occ_csv
      )
    )

  }, ignoreInit = TRUE)

  # --- Réinitialisation ---
  observeEvent(input$reset_btn, {
    taxa$list <- list()
    sugs$df <- NULL
    sug_message("")
    updateTextInput(session, "search", value = "")
    leafletProxy("map") |>
      clearMarkers() |>
      clearPopups()
  })

  # --- Carte ---
  output$map <- renderLeaflet({
    leaflet(
      options = leafletOptions(
        maxZoom = MAX_ZOOM,
        worldCopyJump = TRUE
      )
    ) |>
      setView(lng = 2.5, lat = 46.6, zoom = 6) |>
      addTiles(
        "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
        group = "OpenStreetMap",
        options = tileOptions(maxNativeZoom = MAX_ZOOM)
      ) |>
      addTiles(
        "https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png",
        group = "OpenTopoMap",
        options = tileOptions(maxNativeZoom = 17)
      ) |>
      addTiles(
        "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        group = "Satellite",
        options = tileOptions(maxNativeZoom = 19)
      ) |>
      addLayersControl(
        baseGroups = c("OpenStreetMap", "OpenTopoMap", "Satellite"),
        position = "bottomleft",
        options = layersControlOptions(collapsed = FALSE)
      )
  })

  # --- Affichage des occurrences ---
  observe({
    req(length(taxa$list) > 0)

    proxy <- leafletProxy("map")
    proxy |>
      clearMarkers() |>
      clearPopups()

    for (nm in names(taxa$list)) {
      tx <- taxa$list[[nm]]

      if (is.null(tx$occ) || nrow(tx$occ) == 0) next

      occ <- split(tx$occ, tx$occ$source)

      for (src in names(occ)) {
        d <- occ[[src]]

        shape <- SHAPES[[src]]
        if (is.null(shape)) next

        icon <- leaflet::makeIcon(
          iconUrl = make_icon(shape, tx$color),
          iconWidth = 20,
          iconHeight = 20,
          iconAnchorX = 10,
          iconAnchorY = 10,
          popupAnchorX = 0,
          popupAnchorY = -10
        )

        popup_text <- sprintf(
          paste0(
            "<b>%s</b><br>",
            "Taxon : %s<br>",
            "Source : %s<br>",
            "Observateur : %s<br>",
            "Date : %s<br>",
            "Effectif : %s<br>",
            "Sexe : %s"
          ),
          htmltools::htmlEscape(nm),
          htmltools::htmlEscape(nm),
          htmltools::htmlEscape(if (src == "csv") "CSV PostgreSQL" else if (src == "gbif") "GBIF" else src),
          htmltools::htmlEscape(as.character(d$obs)),
          htmltools::htmlEscape(as.character(d$date)),
          htmltools::htmlEscape(as.character(d$effectif)),
          htmltools::htmlEscape(as.character(d$sexe))
        )

        proxy |>
          addMarkers(
            data = d,
            lng = ~lon,
            lat = ~lat,
            icon = icon,
            group = nm,
            popup = popup_text
          )
      }
    }
  })

  # --- Vignettes latérales ---
  output$taxon_cards <- renderUI({

    if (length(taxa$list) == 0) {
      return(
        tags$p(
          "Aucun taxon sélectionné. Utilisez la barre de recherche.",
          style = "color:#999;"
        )
      )
    }

    lapply(names(taxa$list), function(nm) {
      tx <- taxa$list[[nm]]
      ph <- tx$photos

      img_tag <- if (!is.null(ph) && length(ph) > 0) {
        tags$div(
          class = "photo-gallery",
          tags$div(
            class = "photo-strip",
            lapply(seq_along(ph), function(j) {
              photo <- ph[[j]]
              tags$div(
                class = "photo-slide",
                tags$a(
                  href = photo$url,
                  target = "_blank",
                  rel = "noopener noreferrer",
                  tags$img(
                    src = photo$url,
                    alt = paste("Photographie de", nm),
                    loading = "lazy"
                  )
                ),
                tags$div(
                  class = "credit",
                  HTML(sprintf(
                    "%s &mdash; %s",
                    htmltools::htmlEscape(photo$source),
                    htmltools::htmlEscape(photo$credit)
                  ))
                )
              )
            })
          ),
          if (length(ph) > 1) tags$div(
            class = "credit photo-help",
            sprintf("%d photos — faites défiler horizontalement pour les parcourir", length(ph))
          )
        )
      } else {
        tags$div(
          class = "credit",
          "Aucune photo disponible"
        )
      }

      links <- tags$div(
        class = "taxon-links",
        tags$a(
          href = sprintf(
            "https://www.gbif.org/species/%s",
            tx$usageKey
          ),
          target = "_blank",
          rel = "noopener noreferrer",
          "GBIF"
        ),
        tags$a(
          href = sprintf(
            "https://millibase.myspecies.info/search?search_api_fulltext=%s",
            utils::URLencode(tx$name)
          ),
          target = "_blank",
          rel = "noopener noreferrer",
          "Millibase"
        ),
        tags$a(
          href = sprintf(
            "https://chilobase.biologie.uni-freiburg.de/search?taxon=%s",
            utils::URLencode(tx$name)
          ),
          target = "_blank",
          rel = "noopener noreferrer",
          "ChiloBase"
        )
      )

      n_occ <- if (!is.null(tx$occ)) nrow(tx$occ) else 0

      tags$div(
        class = if (cards_collapsed()) "taxon-card collapsed" else "taxon-card",
        tags$div(
          class = "taxon-title",
          style = sprintf(
            "border-left: 4px solid %s; padding-left: 6px;",
            tx$color
          ),
          nm
        ),
        if (!cards_collapsed()) {
          tags$div(
            class = "taxon-details",
            img_tag,
            links,
            tags$div(
              class = "credit",
              sprintf(
                "%d occurrence(s) affichée(s) — clé GBIF : %s",
                n_occ,
                tx$usageKey
              )
            )
          )
        }
      )
    })
  })
}

shinyApp(ui, server)
