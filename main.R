# Packages and Functiosn --------------------------------------------------
library(tidyverse)
library(rlang)
library(patchwork)
library(stargazer)
library(ivreg)

theme_set(theme_bw())
pal <- RColorBrewer::brewer.pal(8, "Dark2")

prettify_star <- function(star, low = 0.01, high = 1000, digits = 3, scipen = -7) {
  mark  <- '::::'
  
  replace_numbers <- function(x) {
    x <- gsub(mark, '.', x)
    x.num <- as.numeric(x)
    ifelse(
      (x.num >= low) & (x.num < high), 
      round(x.num, digits = digits), 
      prettyNum(x.num, digits = digits, scientific = scipen)
    )
  }    
  
  reg <- paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)")
  gsubfn::gsubfn(reg, ~replace_numbers(x), star)
}

gaze <- function(x, output = "console", ...) {
  star <- stargazer(x,
    covariate.labels = str_replace_all(names(x$coefficients)[-1], 
      c("as\\.factor\\(Siblings\\)" = "Siblings = ")
    ),
    type = if (output == "console") "text" else "latex",
    no.space = TRUE,
    #omit.stat = c("f"),
    df = FALSE,
    omit.table.layout = "n",
    decimal.mark = "::::",
    header = FALSE,
    table.placement = "!htbp",
    ...
  ) %>%
    capture.output()
  
  if (output == "console") {
    cat(star, sep = "\n")
  } else {
    star %>%
      prettify_star() %>%
      writeLines(output)
  }
  
  invisible(star)
}

plot_smooths <- function(data) {
  create_sm <- function(formula, name, data, method = "lm"){
    list(geom_smooth(aes(Siblings, color = !!name), data,
      formula = formula, method = method, se = FALSE, na.rm = TRUE, alpha = 0.7,
      linewidth = 0.75
    ))
  }
  
  list(
    create_sm(y ~ x, "Linear", data = data),
    create_sm(y ~ x + I(x^2), "Quadratic", data = data),
    create_sm(y ~ I(as.factor(trunc(x))), "Dummies", data = data)
  )
}

data_gdensity <- function(...) {
  scalings <- table(data$Siblings) %>% {./max(.)}
  result <- list()
  
  result$data_density <- data %>%
    select(c(Siblings, !!!ensyms(...))) %>%
    group_split(Siblings) %>%
    map_dfr(function(group) {
      scalings_group <- scalings[group$Siblings[1]]
      imap_dfc(select(group, -Siblings), function(var, name) {
        d <- density(var, na.rm = TRUE)
        tibble(
          "{name}==value" := d$x,
          "{name}==dens" := d$y/max(d$y)*0.6*scalings_group
        ) %>%
          arrange("{name}==value")
      }) %>%
        mutate(Siblings = group$Siblings[1], .before = 1)
    }) %>%
    pivot_longer(-Siblings, names_sep = "==", names_to = c("var", ".value"))
  
  result$data_true <- data %>%
    select(c(Siblings, !!!ensyms(...))) %>%
    mutate(Siblings_jitter = Siblings + runif(n(), 0, 0.6)) %>%
    pivot_longer(-starts_with("Siblings"), names_to = c("var"))
  
  result
}

plot_gdensity <- function(graph_data, coord, geom = "poly") {
  other_coord <- c(x = "y", y = "x")[coord]
  aes_list <- list(
    exprs("{coord}" := Siblings_jitter),
    exprs("{coord}" := Siblings + dens, group = Siblings)
  )
  list(
    geom_point(aes(!!!aes_list[[1]]), graph_data$data_true, alpha = 0.05),
    if (geom == "poly") {
      geom_polygon(aes(!!!aes_list[[2]]),
        linewidth = 0.5, color = pal[7], fill = pal[7], alpha = 0.5
      )
    } else {
      geom_path(aes(!!!aes_list[[2]]), linewidth = 0.5, color = pal[7])
    }
  )
}

mean_wona <- function(x) {
  x %>% replace_na(0) %>% mean()
}

quick_agg <- function(..., filter = TRUE) {
  columns <- exprs(...)
  
  result <- data_raw %>%
    group_by(SERIAL) %>%
    summarise(
      Siblings = sum(!!is_child, na.rm = TRUE),
      !!!columns
    ) %>%
    {if (filter) filter(., Siblings > 0 & Siblings < 7) else .} %>%
    select(-(1:2))
  
  print(map(result, ~ round(table(.x)/length(.x), 4) * 100))
  invisible(result)
}

plot_partial <- function(updt, iv) {
  data_partial <- map_dfr(model3$formulas[1:2], ~ tibble(
    Siblings = lm(update(.x, updt), data)$residuals,
    y = lm(update(update(.x, Siblings ~ .), updt), data)$residuals,
    model = as.character(.x[[2]])
  ))
  
  ggplot(data_partial, aes(Siblings, y)) +
    geom_point(alpha = 0.05) +
    geom_smooth(method = "lm", color = pal[1]) +
    facet_wrap(vars(model), scales = "free") +
    labs(
      title = "Partial Regression of Quality on Siblings",
      subtitle = c("With Controls", "With Instruments")[iv + 1],
      y = "Residuals of Y", x = "Residuals of Siblings"
    )
}



# Data --------------------------------------------------------------------
# Importing and removing fixed variables; Saving as RDS:
if (FALSE) {
  data_raw <- read.csv("Data/ipumsi_data.csv") %>%
    as_tibble()
  
  data_raw <- data_raw %>%
    select(-all_of(map_dbl(data_raw, ~length(unique(.x))) %>% {names(.)[. == 1]}))
  
  saveRDS(data_raw, "Data/data_raw.RDS")
}

# Filtering out negative incomes; transforming missing codes to NA, etc.:
data_raw <- readRDS("Data/data_raw.RDS") %>%
  filter(INCTOT >= 0) %>%
  mutate(
    across(where(is.numeric) & !SERIAL,
      ~ifelse(grepl("99+8?", as.character(.x)), NA, .x)
    ),
    STEPPOP = c("0" = 0, "2" = 3, "3" = 5)[as.character(STEPPOP)],
    RACE = ifelse(RACE %in% 40:49, 40, RACE)
  )

# Setup:
is_child <- expr(RELATE == 3 & (AGE <= 17 | (AGE <= 20 & MARST == 1 & is.na(HRSUSUAL1))))
is_parent <- expr(RELATE %in% 1:2)

# Aggregating data to the household level. Creating relevant variables:
data <- data_raw %>%
  filter(NFAMS == 1 & GQ == 10) %>%
  arrange(AGE) %>%
  group_by(HouseID = as.factor(SERIAL)) %>%
  summarise(
    #-- Technical -¬
    PresentMom = all(MOMLOC[!!is_child] != 0),
    PresentPop = all(POPLOC[!!is_child] != 0),
    Siblings = sum(!!is_child, na.rm = TRUE),
    ChildAges = list(AGE[!!is_child]),
    ChildMales = list(SEX[!!is_child] == 1),
    PernumMom = list(PERNUM_MOM[!!is_child]),
    HasOldChild = any(AGE[RELATE == 3] >= 21),
    HasGrandchilds = any(RELATE == 4, na.rm = TRUE),
    #-- Quality measures -¬
    BedroomsCapta = BEDROOMS[1]/PERSONS[1],
    YearsEduc = mean_wona(AGE[!!is_child] - 5 - YRSCHOOL[!!is_child]),
    StudentsCapta = sum(SCHOOL[!!is_child] == 1)/Siblings,
    HrsUsual = list(HRSUSUAL1[!!is_child] %>% ifelse(is.na(.), 0, .)),
    #-- Controls -¬
    MaleAvg = mean_wona(SEX[!!is_child] == 1),
    HasDisabled = (DISABLED[!!is_child][1] == 2) && any(DISABLED[!!is_child] == 1, na.rm = TRUE),
    FirstDisabled = DISABLED[!!is_child] %>% {.[length(.) - 1] == 1},
    FamSize = sum(!RELATE %in% 1:4, na.rm = TRUE),
    HasStep = any(STEPMOM[!!is_child] == 3) | any(STEPPOP[!!is_child] == 3),
    ParentsAgeAvg = mean_wona(c(AGE[!!is_parent], YRSCHOOL_MOM, YRSCHOOL_POP)),
    ParentsEducAvg = mean_wona(c(YRSCHOOL[!!is_parent], YRSCHOOL_MOM[!!is_child], YRSCHOOL_POP[!!is_child])),
    ParentsWorkAvg = mean_wona(c(HRSUSUAL1[!!is_parent], HRSUSUAL1_MOM[!!is_child], HRSUSUAL1_POP[!!is_child])),
    ParentsRace = RACE[!!is_parent][1] %>%
      {ifelse((sum(!!is_parent) > 1) && (. != RACE[!!is_parent][2]), 0, .)} %>% as.factor(),
    ParentsUnited = MARST[!!is_parent][1] == 2, #& MARST_MOM == 2 & MARST_POP == 2
    ParentsCitizenship = CITIZEN[!!is_parent][1] %>%
      {ifelse((sum(!!is_parent) > 1) && (. != CITIZEN[!!is_parent][2]), 0, .)} %>% as.factor(),
    #-- Income -¬
    IncTot = sum(INCTOT, na.rm = TRUE),
    IncWel = if (IncTot > 0) sum(INCWEL, na.rm = TRUE)/IncTot else 0,
    IncPen = if (IncTot > 0) sum(INCRET, na.rm = TRUE)/IncTot else 0,
  ) %>%
  ungroup() %>%
  filter(Siblings > 0 & Siblings < 7 & PresentMom & !HasOldChild & !HasGrandchilds) %>%
  mutate(
    IncTotQuant = ecdf(IncTot)(IncTot),
    IncTotLog = log(IncTot + 1, base = 10),
    IncTotCut = cut(IncTot, c(-Inf,25000,50000,100000,Inf)) %>%
      fct_relabel(., \(x) set_names(paste(1:4), levels(.))[x])
  ) %>%
  rowwise() %>%
  mutate(
    HasWorkers = any(HrsUsual != 0),
    WorkedCapta = sum(HrsUsual)/Siblings,
    ChildsAgeMin = min(ChildAges, na.rm = TRUE), #AgeMax = max(ChildAges, na.rm = TRUE),
    ChildsAgeAvg = mean_wona(ChildAges),
    HasTwins = any(map2_lgl(duplicated(ChildAges), duplicated(PernumMom), `&`)) &
      duplicated(ChildAges)[length(ChildAges)],
    OldTwins = HasTwins & ChildAges %>% {.[which.min(.)] >= 7},
    FirstSameSex = map2_lgl(duplicated(ChildMales[1:2]), duplicated(PernumMom)[1:2], `&`)[2],
    TwinsFirst = HasTwins & length(ChildAges) == 2,
    TwinsSecond = HasTwins & length(ChildAges) == 3,
    TwinsThird = HasTwins & length(ChildAges) == 4,
  ) %>%
  ungroup() %>%
  na.omit() #check NA's with: map_dbl(data, ~sum(is.na(.x)))



# Exploratory Analysis ----------------------------------------------------
# -------------------- Filtering --------------------
# Analyzing NA and NaN values:
data %>%
  select(-c(ChildAges, ChildMales, HrsUsual, PernumMom)) %>%
  map_df(~ c(sum(is.nan(.x)), sum(are_na(.x)))) %>%
  print(width = 9999)


# Households and households with children
data_raw$SERIAL %>% unique() %>% length()
filter(data_raw, RELATE == 3)$SERIAL %>% unique() %>% length()

quick_agg(
  Under17 = sum(RELATE == 3 & AGE <= 17, na.rm = TRUE) == 0,
  Under21 = sum(RELATE == 3 & (AGE <= 17 | (AGE <= 20 & MARST == 1 & is.na(HRSUSUAL1))), na.rm = TRUE) == 0,
  filter = FALSE
) %>%
  pull(Under21) %>%
  table()


# Families per household:
quick_agg(NFams = NFAMS[1] > 1, Gq = GQ[1] != 10)


# Location and step parenting:
quick_agg(
  MomLoc = any(MOMLOC[!!is_child] != 0),
  PopLoc = any(POPLOC[!!is_child] != 0),
  Both = any((MOMLOC[!!is_child] != 0) & (POPLOC[!!is_child] != 0)),
  Any = any((MOMLOC[!!is_child] != 0) | (POPLOC[!!is_child] != 0))
)

quick_agg(Step = any((STEPMOM[!!is_child] != 0) | (STEPPOP[!!is_child] != 0)))


# Old children and grandchildren:
quick_agg(Over21 = any(AGE[RELATE == 3] >= 21))

data_raw %>%
  filter(RELATE == 3) %>%
  pull(AGE) %>%
  ggpubr::gghistogram(ggtheme = theme_bw()) +
  labs(title = "Histogram of 'Childs' Ages", y = "Count", x = "Age")

ggsave("Figures/ages.png", width = 16, height = 12, units = "cm")

quick_agg(HasGrandchilds = any(RELATE == 4, na.rm = TRUE))

# -------------------- Controls --------------------
# Twins:
data %>%
  select(contains("Twins")) %>%
  map(table)


# Disabled:
data$HasDisabled %>% table()


# Years of education in children:
data_raw %>%
  filter(AGE < 21) %>%
  select(AGE, YRSCHOOL, EDUCPR) %>%
  pivot_longer(-AGE) %>%
  ggplot(aes(AGE, value)) +
  geom_count() +
  geom_abline(aes(alpha = name, slope = 1, intercept = -5), color = pal[1], linewidth = 1) +
  scale_alpha_manual(values = c(0, 1)) +
  theme(legend.position = "none") +
  facet_wrap(vars(name), scales = "free_y") +
  labs(title = "Education of Children by Age", y = "Value", x = "Age")

ggsave("Figures/educ.png", width = 16, height = 12, units = "cm")

# Analizyng old children:
data_raw %>%
  filter(RELATE == 3 & AGE <= 20) %>%
  pull(MARST) %>%
  {round(table(.)/length(.), 4) * 100}

data %>%
  rowwise() %>%
  mutate(
    HasWorkersSub18 = map2_lgl(HrsUsual[[1]] != 0, ChildAges[[1]] < 18, `&`),
    HasWorkersSub16 = map2_lgl(HrsUsual[[1]] != 0, ChildAges[[1]] < 16, `&`)
  ) %>%
  select(starts_with("HasWorkers")) %>%
  map(table)

map(list("16 to 17" = 16:17, "18 to 20" = 18:20), function(ages) {
  data_raw %>%
    filter(RELATE == 3 & AGE %in% ages) %>%
    pull(EDUCPR) %>%
    {c(sum(. < 410), sum(. == 410), sum(. > 410))/length(.)} %>%
    round(3)
})



# Task 1 ------------------------------------------------------------------
# -------------------- Models --------------------
model1 <- list()

model1$formulas <- list(
  "Linear" = Siblings ~ IncTot,
  "Linear + bins" = Siblings ~ IncTot + IncTotCut,
  "Linear * bins" = Siblings ~ IncTot*IncTotCut,
  "Quadratic" = Siblings ~ IncTot + I(IncTot^2),
  "Log" = Siblings ~ IncTot + IncTotLog,
  "Quantile" = Siblings ~ IncTotQuant,
  "Quantile * bins" = Siblings ~ IncTotQuant*IncTotCut,
  "Quantile Q." = Siblings ~ IncTotQuant + I(IncTotQuant^2)
)

model1$models <- map(model1$formulas, ~lm(.x, data))

model1$results <- imap_dfr(model1$models, function(m, name) {
  tibble(model = name, fit = fitted(m), res = residuals(m)) %>%
    mutate(
      value = if (name %in% c("Quantile", "Quantile * bins", "Quantile Q.")) {
        quantile(data$IncTot, m$model[,"IncTotQuant"])
      } else if (name %in% c("Log")) {
        10^m$model[,"IncTotLog"] - 1
      } else {
        m$model[,"IncTot"]
      },
      model_cat = `if`(grepl("Linear", name), "Linear", `if`(
        grepl("Quantile", name), "Quantile", "Non-Linear"
      )),
      model = factor(model)
    )
})

model1$models %>%
  .[grepl("Linear", names(.))] %>%
  gaze(title = "Linear CEFs", output = "Tables/cefs_linear.tex")

model1$models %>%
  .[grepl("Quadratic|Log", names(.))] %>%
  gaze(title = "Non-linear CEFs", output = "Tables/cefs_non-linear.tex")

model1$models %>%
  .[grepl("Quantile", names(.))] %>%
  gaze(title = "Quantile CEFs", output = "Tables/cefs_quantile.tex")

# -------------------- Graphs --------------------
graph1 <- data_gdensity(IncTot)

ggplot(graph1$data_density, aes(value, Siblings)) +
  plot_gdensity(graph1, "y") +
  geom_line(aes(y = fit, color = model), model1$results,
    linewidth = 0.75, alpha = 0.7
  ) +
  facet_wrap(vars(model_cat), nrow = 1) +
  scale_color_manual(values = c(pal[1:6], "darkred", "darkblue"), drop = FALSE) +
  ylim(2,5.9) +
  labs(
    title = "Siblings versus Total Income",
    x = "Total Income",
    y = "Siblings",
    color = "Model"
  ) +
  theme(legend.position = "bottom")

ggsave("Figures/cefs.png", width = 16, height = 12, units = "cm")


# Task 2 ------------------------------------------------------------------
# -------------------- Models --------------------
model2 <- list()

model2$formulas <- expand.grid(
  c("BedroomsCapta", "YearsEduc", "StudentsCapta", "HasWorkers", "WorkedCapta"),
  "~",
  c("Siblings", "as.factor(Siblings)")
) %>%
  apply(1, \(x) as.formula(paste(x, collapse = " ")))

model2$models <- map(model2$formulas, ~lm(.x, data))

model2$models[c(1, 2, 6, 7)] %>%
  gaze(title = "Siblings Effect - no Controls", output = "Tables/quality.tex")


# -------------------- Graphs --------------------
graph2 <- data_gdensity(BedroomsCapta, YearsEduc, StudentsCapta, WorkedCapta)

graph2$data_hasworkers <- data %>%
  select(c(Siblings, HasWorkers)) %>%
  mutate(HasWorkers = as.integer(HasWorkers), var = "HasWorkers")

graph2$graph_cont <- ggplot(graph2$data_density, aes(Siblings, value)) +
  plot_gdensity(graph2, "x") +
  plot_smooths(graph2$data_true) +
  facet_wrap(vars(var), scales = "free_y") +
  scale_color_manual(values = pal) +
  labs(y = "Value", color = "Models")

graph2$graph_disc <- ggplot(graph2$data_hasworkers, aes(Siblings, HasWorkers)) +
  geom_count(color = pal[7]) +
  plot_smooths(graph2$data_hasworkers) +
  facet_wrap(vars(var)) +
  scale_size_continuous(breaks = c(100, 1000)) +
  scale_color_manual(values = pal) +
  labs(x = "Siblings", y = "", size = "Count/Density", color = "Models")

graph2$layout <- "AAAAB\nAAAAB\nAAAAC\nAAAAC"

graph2$graph_cont + graph2$graph_disc + guide_area() +
  plot_layout(widths = c(4,1), guides = "collect", design = graph2$layout) +
  plot_annotation(title = "Quality Measures versus Siblings")

ggsave("Figures/quality.png", width = 16, height = 14, units = "cm")


# Task 3 ------------------------------------------------------------------
# -------------------- Models --------------------
model3 <- list()

model_controls <- " + MaleAvg + ChildsAgeAvg + ChildsAgeMin + HasDisabled +
     FamSize + IncTotQuant + IncWel + IncPen +
     ParentsAgeAvg + ParentsEducAvg + ParentsWorkAvg + ParentsUnited +
     ParentsRace + ParentsCitizenship + HasStep" %>%
  str_remove_all("\n   ")

model3$formulas <- expand.grid(
  quality = c("BedroomsCapta", "YearsEduc"),
  siblings = c("~ Siblings", "~ as.factor(Siblings)"),
  controls = model_controls
) %>%
  apply(1, \(x) as.formula(paste(x, collapse = " ")))

model3$models <- map(model3$formulas, ~ lm(.x, data))
map(model3$models, car::vif)

gaze(model3$models,
     keep = "Siblings",
     title = "Siblings Effect - with Controls",
     output = "Tables/quality_controls.tex"
)


# -------------------- Graphs --------------------
graph3 <- list()

graph3$data_cor <- data %>%
  select(c( #ParentsRace, ParentsCitizenship
    MaleAvg, ChildsAgeAvg, ChildsAgeMin, HasDisabled,
    FamSize, IncTotQuant, IncWel, IncPen, 
    ParentsAgeAvg,  ParentsEducAvg, ParentsWorkAvg, ParentsUnited, HasStep
  )) %>% 
  cor() %>%
  as_tibble(rownames = "var1") %>%
  pivot_longer(-var1, names_to = "var2")

ggplot(graph3$data_cor, aes(var1, var2, fill = value)) +
  geom_tile() +
  labs(
    title = "Correlation of Controls", x = "Variable",  y = "Variable",
    fill = "Correlation"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8))
ggsave("Figures/correlation.png", width = 16, height = 12, units = "cm")

plot_partial(~ . - Siblings, FALSE)
ggsave("Figures/quality_controls.png", width = 16, height = 12, units = "cm")


# Task 4 ------------------------------------------------------------------
# -------------------- Models --------------------
model4 <- list()

model4$formulas <- expand.grid(
  quality = c("BedroomsCapta", "YearsEduc"),
  basic1 = paste("~ Siblings", model_controls),
  basic2 = "|",
  instruments = c("HasTwins", "TwinsSecond", "TwinsThird", "FirstSameSex", "FirstDisabled"),
  controls = model_controls
) %>%
  arrange(quality) %>%
  apply(1, \(x) as.formula(paste(x, collapse = " ")))

model4$models <- map(model4$formulas, ~ ivreg(.x, data = data))

gaze(model4$models[1:5], output = "Tables/rooms_instrument.tex",
     keep = "Siblings",
     title = "Siblings Effect on Bedrooms - with Instrument"
)

gaze(model4$models[6:10], output = "Tables/school_instrument.tex",
  keep = "Siblings",
  title = "Siblings Effect on Schooling - with Instrument"
)


# -------------------- Graphs --------------------
graph4 <- list()

map_dbl(
  select(data, c(HasTwins, FirstSameSex, FirstDisabled, TwinsSecond, TwinsThird)),
  ~ cor(.x, data$Siblings) %>% round(3) %>% `*`(100) 
)

#plot_partial(~ . - Siblings + TwinsThird, TRUE)
#ggsave("Figures/quality_instruments.png", width = 16, height = 12, units = "cm")


# Other tasks -------------------------------------------------------------
# Correlations:
map_dbl(
  c("Siblings", "MaleAvg", "ChildsAgeAvg", "ChildsAgeMin", "HasDisabled",
    "FamSize", "IncTotQuant", "IncWel", "IncPen", "ParentsAgeAvg", "ParentsEducAvg",
    "ParentsWorkAvg", "ParentsUnited", "HasStep") %>% set_names(),
  ~ round(cor(data$IncTot, data[[.x]]), 3)
)

map_dbl(
  c("IncTot", "IncTotQuant", "MaleAvg", "ChildsAgeAvg", "ChildsAgeMin", "HasDisabled",
    "FamSize", "IncTotQuant", "IncWel", "IncPen", "ParentsAgeAvg", "ParentsEducAvg",
    "ParentsWorkAvg", "ParentsUnited", "HasStep") %>% set_names(),
  ~ round(cor(data$Siblings, data[[.x]]), 3)
)

# Diagnostics:
walk(
  list(
    model1$models$`Quantile Q.`,
    model2$models[[2]],
    model3$models[[2]],
    model4$models[[2]]
  ),
  ~ lmtest::bptest(.x) %>% broom::tidy() %>% stargazer(summary = FALSE)
)

#coeftest(food.ols, vcov = vcovHC(food.ols, "HC1"))
