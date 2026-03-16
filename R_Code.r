library(tidyverse)
library(maps)
library(sf)
library(usmap)
library(gganimate)
library(gifski)
library(tigris)
library(prism)
prism_set_dl_dir("~/prism")
library(exactextractr)
library(raster)
library(terra)
library(stringr)
library(pROC)
library(PRROC)
library(xgboost)
library(zoo)
library(randomForest)

wildfires <- read.csv("C:/Users/krimi/Downloads/wildfires_sample_100k.csv")
str(wildfires)

#we want to find a binary variable (aka either "yes" or "no") that can indicate if a notable wildfire occurs (one spatial one temporal). 
# To discover time patterns, we can create FIRE_DATE to better rep time
wildfires <- wildfires %>% filter(!(is.na(FIPS_CODE)))
wildfires <- wildfires %>% mutate(FIRE_DATE = as.Date(DISCOVERY_DOY-1, origin = paste0(FIRE_YEAR, "-01-01"))) %>% mutate(MONTH = month(FIRE_DATE)) %>% mutate(MONTH = factor(MONTH, levels=1:12, labels = month.name)) %>% mutate(statefips = fips(STATE), FIPS_CODE = sprintf("%03d", FIPS_CODE))

#we'll filter and cut off at recommended cutoff firesize >= 300
notable_wildfires <- wildfires %>% filter(FIRE_SIZE >= 300)

#we'll organize by counties and fire size to analyze the top 5 years with the largest wildfires and where they were
counties_burned <- wildfires %>% group_by(statefips, FIPS_CODE, FIRE_YEAR) %>% summarize(acres_burned = sum(FIRE_SIZE, na.rm = TRUE), .groups="drop") 
counties_sf <- counties(cb = TRUE, class = "sf")
county_map <- counties_sf %>% left_join(counties_burned, by = c("COUNTYFP" = "FIPS_CODE", "STATEFP" = "statefips"))

top_5_years <- counties_burned %>%
  group_by(FIRE_YEAR) %>%
  summarize(total_per_year = sum(acres_burned, na.rm = TRUE)) %>%
  arrange(desc(total_per_year)) %>%
  slice_head(n = 5) %>%
  pull(FIRE_YEAR)
county_map_filtered <- county_map %>% filter(!is.na(FIRE_YEAR)) %>% filter(FIRE_YEAR %in% top_5_years)

county_map_plot <- ggplot() +
  geom_sf(data = counties_sf, fill = "white", color = "darkgray", size = 0.05) +
  geom_sf(data = county_map_filtered, aes(fill = acres_burned), color = NA) +
  scale_fill_viridis_c(trans = "log10", na.value = "transparent") +
  coord_sf(xlim = c(-125, -65), ylim = c(24, 50), expand = FALSE) +
  facet_wrap(~ FIRE_YEAR, ncol=2) +
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5), legend.title = element_text(size = 25, face = "bold"), legend.text = element_text(size = 18), legend.key.width = unit(2, "cm"), strip.text = element_text(size=25, face= "bold")) + 
  labs(title = "Acres Burned by County in Top 5 Years with Largest Wildfires", fill = "Acres Burned (Log10)")


#For this plot, there seems to be a tendency for large wildfires to happen on the west side of the US during the late 2000s and early 2010s.

# For second plot, I will explore how seasons affect where wildfires are. We'll select the year with the most wildfires and analyze how each passing month affects where the wildfires are located.

notable_wildfires %>% group_by(FIRE_YEAR) %>% summarize(top_fire = n()) %>% arrange(desc(top_fire)) %>% head(1)
#2006 had the most wildfires
most_fires <- notable_wildfires %>% filter(FIRE_YEAR == 2006)
fires_monthly <- most_fires %>% group_by(MONTH, COUNTY) %>% summarize(acres_burned = sum(FIRE_SIZE), .groups = "drop")

monthly_map <- county_map %>% left_join(fires_monthly, by = c("NAME" = "COUNTY"))
monthly_map_filtered <- monthly_map %>% filter(FIRE_YEAR == 2006, !is.na(MONTH)) 

monthly_fires_plot <- ggplot() + geom_sf(data=counties_sf, fill = "white", color = "gray", size = 0.05) + geom_sf(data=monthly_map_filtered, aes(fill=acres_burned.x), color = "black") + scale_fill_viridis_c() + coord_sf(xlim = c(-125, -65), ylim = c(24, 50), expand = FALSE) + labs(title = "Plot A: Notable Wildfires in 2006", subtitle = "Month: {closest_state}", fill = "Acres Burned") + transition_states(MONTH, transition_length=2, state_length=1) 

anim <- gganimate::animate(monthly_fires_plot, nframes=12, fps=2, renderer = gifski_renderer())
anim_save("wildfire_animation.gif", animation=anim)
#animation tells us that in the spring and summer, wildfires tend to be on the west side while in fall and winter, they tend to be in the midwest/east

#we'll construct the binary variable as 1 if a fire thats size was > 300 occured in that month and county and 0 if not
model_data <- wildfires %>% filter(!(is.na(FIRE_SIZE))) %>% mutate(notable_fire = ifelse(FIRE_SIZE >= 300, 1, 0)) %>% dplyr::select(FIRE_YEAR, FIRE_SIZE, STATE, FIPS_CODE, MONTH, FIRE_DATE, LATITUDE, LONGITUDE, notable_fire) %>% filter(!(is.na(FIPS_CODE)))

# attatch prism variables to data
types = c("ppt", "tmean")
years <- 1992:2015
months <- 1:12
for (t in types){
  get_prism_monthlys(type = t, years = years, mon = months, keepZip = FALSE)
}
climate_stack <- pd_stack(prism_archive_subset("tmean", "monthly", resolution="4km"))
climate_stack2 <- pd_stack(prism_archive_subset("ppt", "monthly", resolution="4km"))
climate_rast <- rast(climate_stack)
climate_rast2 <- rast(climate_stack2)

wildfire_sf <- st_as_sf(model_data, coords = c("LONGITUDE", "LATITUDE"), crs = 4269)
wildfire_sf <- st_transform(wildfire_sf, crs(climate_rast))
extracted_climate <- terra::extract(climate_rast, wildfire_sf)
wildfire_sf <- st_transform(wildfire_sf, crs(climate_rast2))
extracted_ppt <- terra::extract(climate_rast2, wildfire_sf)
wildfires_with_climate <- bind_cols(wildfire_sf, extracted_climate, extracted_ppt)

# formatting properly
wildfires_final <- wildfires_with_climate %>%
  pivot_longer(
    cols = starts_with("prism_"), 
    names_to = "climate_var", 
    values_to = "value"
  ) %>%
  # extract the last 6 digits (the YYYYMM part) from the column name
  mutate(
    var_type = if_else(str_detect(climate_var, "tmean"), "temp", "ppt"),
    date_code_prism = str_extract(climate_var, "\\d{6}$"),
    # create a matching code from fire data (ex: 1998 and 4 becomes "199804")
    date_code_fire = sprintf("%d%02d", FIRE_YEAR, MONTH)
  ) %>% filter(date_code_prism == date_code_fire) %>%
  dplyr::select(-climate_var) %>% ungroup()
wildfires_final <- wildfires_final %>% pivot_wider(names_from = var_type, values_from = "value") 
wildfires_final <- wildfires_final %>%
  mutate(
    FIPS_CODE = sprintf("%03d", as.numeric(FIPS_CODE)),  # ensure 3-digit string
    FIPS_CODE = factor(FIPS_CODE),
    STATE     = factor(STATE)
  )

#cleaning up and creating variables to prevent data leakage
wildfires_final <- wildfires_final %>%
  arrange(FIPS_CODE, FIRE_YEAR, MONTH) %>%
  group_by(FIPS_CODE) %>%
  mutate(ppt_3mo = rollsum(ppt, 3, fill = NA, align = "right"), temp_3mo = rollmean(temp, 3, fill = NA, align = "right"), FIRE_SIZE_3mo = rollmean(FIRE_SIZE, 3, fill = NA, align = "right"), fire_lag12  = lag(FIRE_SIZE, 12), fire_roll3  = rollsum(FIRE_SIZE, 3, fill = NA, align = "right"), fire_roll12 = rollsum(FIRE_SIZE, 12, fill = NA, align = "right"), size_roll3  = rollmean(FIRE_SIZE, 3, fill = NA, align = "right"), size_roll12 = rollmean(FIRE_SIZE, 12, fill = NA, align = "right"), y_next = lead(notable_fire, 1)) %>% ungroup() %>%
filter(!is.na(ppt_3mo), !is.na(temp_3mo), !is.na(FIRE_SIZE_3mo), !is.na(y_next), !is.na(temp),!is.na(ppt),!is.na(MONTH),!is.na(FIPS_CODE),!is.na(STATE),!is.na(FIRE_SIZE),!is.na(FIRE_YEAR),!is.na(notable_fire))


wildfires_final <- wildfires_final %>%
  mutate(ppt = as.numeric(ppt), temp = as.numeric(temp), ppt_3mo = as.numeric(ppt_3mo), temp_3mo = as.numeric(temp_3mo), FIRE_SIZE_3mo = as.numeric(FIRE_SIZE_3mo), FIPS_CODE = factor(FIPS_CODE), STATE = factor(STATE), MONTH = factor(MONTH),
y_next = factor(y_next, levels = c(0,1)))

coords <- st_coordinates(wildfires_final)
wildfires_final$LON <- coords[,1]
wildfires_final$LAT <- coords[,2]
wildfires_final <- st_drop_geometry(wildfires_final)
years_sorted <- sort(unique(wildfires_final$FIRE_YEAR))

#mini table to analyze who has more acres burned
wildfires_plot <- wildfires_final %>%
  mutate(region = ifelse(LON < -100, "West", "East"))

wildfires_plot <- wildfires_plot %>%
  group_by(region) %>%
  summarise(total_acres = sum(FIRE_SIZE, na.rm = TRUE))
saveRDS(wildfires_plot, file = "region_acres_burned.rds")

# Summary table of all predictors and summary stats
wildfire_summary <- wildfires_final %>%
  summarise(
    n_rows = n(),
    n_fires = sum(notable_fire, na.rm = TRUE),
    spatial_unit = "FIPS_CODE",
    temporal_unit = "MONTH",
    predictors = paste(c("temp_3mo", "ppt_3mo", "LAT", "LON", "STATE", "FIRE_SIZE_3mo", "FIRE_YEAR", "FIPS_CODE", "MONTH"), 
    collapse = ", "), mean_temp = mean(temp, na.rm = TRUE), sd_temp = sd(temp, na.rm = TRUE), mean_ppt = mean(ppt, na.rm = TRUE),
    sd_ppt = sd(ppt, na.rm = TRUE), mean_fire_size = mean(FIRE_SIZE, na.rm = TRUE), sd_fire_size = sd(FIRE_SIZE, na.rm = TRUE),        
    prop_notable = mean(notable_fire, na.rm = TRUE), prop_next_month = mean(as.numeric(y_next) - 1, na.rm = TRUE))


# plot of seasonality in the year with the most wildfires
top_year <- wildfires_final %>%
  group_by(FIRE_YEAR) %>%
  summarize(total_y = sum(as.numeric(as.character(y_next))), .groups = "drop") %>%
  arrange(desc(total_y)) %>%
  slice(1) %>%
  pull(FIRE_YEAR)

top_year_data <- wildfires_final %>%
  filter(FIRE_YEAR == top_year)

monthly_y <- top_year_data %>%
  group_by(MONTH) %>%
  summarize(y_next_total = sum(as.numeric(as.character(y_next))), y_next_mean  = mean(as.numeric(as.character(y_next))),
            .groups = "drop") 
monthly_y <- monthly_y %>% mutate(MONTH_label = factor(MONTH, levels = month.name))

y_monthly_plot <- ggplot(monthly_y, aes(x = MONTH_label, y = y_next_mean, group = 1)) +
  geom_line(color = "firebrick", size = 1.2) +
  geom_point(color = "firebrick", size = 3) +
  labs(
    title = paste("Plot B: Seasonality of y_next in", top_year), subtitle = "The year with the most wildfires",
    x = "Month",
    y = "Average y_next (probability of notable fire)"
  ) + transition_reveal(as.numeric(MONTH_label)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

anim <- gganimate::animate(y_monthly_plot, nframes = 12, fps = 2, renderer = gifski_renderer())
anim_save("y_next_monthly_animation.gif", animation=anim)

#formula all models will use
frmla <- y_next ~ temp_3mo + ppt_3mo + MONTH + LAT + LON + STATE + FIRE_SIZE_3mo + fire_lag12 + fire_roll3 + fire_roll12 + size_roll3 + size_roll12
min_train_years <- 3
folds <- tibble(
  train_end_year = years_sorted[(min_train_years):(length(years_sorted) - 2)],
  valid_year = years_sorted[(min_train_years + 1):(length(years_sorted) - 1)]
)

folds <- folds %>%
  mutate(
    # GLM MODEL
    auc_glm = map2_dbl(train_end_year, valid_year, function(tr_end, va_year) {

      train <- wildfires_final %>% 
        filter(FIRE_YEAR <= tr_end) %>%
        filter(!is.na(FIPS_CODE), !is.na(STATE), !is.na(MONTH))

      valid <- wildfires_final %>% 
        filter(FIRE_YEAR == va_year) %>%
        filter(as.character(FIPS_CODE) %in% as.character(train$FIPS_CODE),
               STATE %in% train$STATE) %>%
        filter(!is.na(FIPS_CODE), !is.na(STATE), !is.na(MONTH))

      fit <- glm(frmla, data = train, family = binomial)
      p_hat <- predict(fit, newdata = valid, type = "response")
      roc_obj <- roc(response = valid$y_next, predictor = p_hat, quiet = TRUE)
      as.numeric(auc(roc_obj))}),

    # RANDOM FOREST MODEL
    auc_rf = map2_dbl(train_end_year, valid_year, function(tr_end, va_year) {
      wildfires_rf <- wildfires_final %>%
  mutate(fire_lag12 = replace_na(fire_lag12, 0), fire_roll3 = replace_na(fire_roll3, 0), fire_roll12 = replace_na(fire_roll12, 0),
    size_roll3 = replace_na(size_roll3, 0), size_roll12 = replace_na(size_roll12, 0))
      train <- wildfires_rf %>% 
        filter(FIRE_YEAR <= tr_end) %>%
        filter(!is.na(FIPS_CODE), !is.na(STATE), !is.na(MONTH))

      valid <- wildfires_rf %>% 
        filter(FIRE_YEAR == va_year) %>%
        filter(as.character(FIPS_CODE) %in% as.character(train$FIPS_CODE),
               STATE %in% train$STATE) %>%
        filter(!is.na(FIPS_CODE), !is.na(STATE), !is.na(MONTH)) %>% drop_na()
      rf_fit <- randomForest(frmla, data = train)
      p_hat <- predict(rf_fit, newdata = valid, type = "prob")[,2]
      roc_obj <- roc(response = valid$y_next, predictor = p_hat, quiet = TRUE)
      as.numeric(auc(roc_obj))}))

folds %>%
  summarise(
    mean_glm_auc = mean(auc_glm, na.rm = TRUE),
    mean_rf_auc = mean(auc_rf, na.rm = TRUE)
  )

# formula and features
features <- c("temp_3mo","ppt_3mo","FIRE_SIZE_3mo","LAT","LON","STATE","fire_lag12", "fire_roll3", "fire_roll12", "size_roll3","size_roll12","MONTH", "FIRE_YEAR", "FIPS_CODE")
target <- "y_next"

# Ensure factors are converted to integer codes consistently
encode_features <- function(df, train_levels = NULL) {
  df <- df %>% mutate(across(c("STATE","FIPS_CODE","MONTH"), as.factor))
  if(!is.null(train_levels)){
    for(col in c("STATE","FIPS_CODE","MONTH")){
      df[[col]] <- factor(df[[col]], levels = train_levels[[col]])
      df[[col]] <- as.integer(df[[col]])
    }
  } else {
    for(col in c("STATE","FIPS_CODE","MONTH")){
      df[[col]] <- as.integer(df[[col]])
    }
  }
  df
}

folds2 <- tibble(
  train_end_year = years_sorted[(min_train_years):(length(years_sorted) - 2)],
  valid_year     = years_sorted[(min_train_years + 1):(length(years_sorted) - 1)]
)

folds2 <- folds2 %>%
  mutate(
    auc_xgb = map2_dbl(train_end_year, valid_year, function(tr_end, va_year) {
      
      # Training and validation splits (time-aware)
      train <- wildfires_final %>% filter(FIRE_YEAR <= tr_end) %>% drop_na(all_of(features), all_of(target))
      valid <- wildfires_final %>% filter(FIRE_YEAR == va_year) %>% drop_na(all_of(features), all_of(target))
      
      # Save factor levels from training
      train_levels <- lapply(train[,c("STATE","FIPS_CODE","MONTH")], levels)
      
      # Encode categorical features as integers
      train <- encode_features(train)
      valid <- encode_features(valid, train_levels)
      
      # Matrix for XGBoost
      X_train <- as.matrix(train[,features])
      y_train <- as.numeric(train[[target]]) - 1
      X_valid <- as.matrix(valid[,features])
      y_valid <- as.numeric(valid[[target]]) - 1
      
      # Compute scale_pos_weight
      scale_pos_weight <- sum(y_train == 0) / sum(y_train == 1)
      
      # Create DMatrix
      dtrain <- xgb.DMatrix(data = X_train, label = y_train)
      dvalid <- xgb.DMatrix(data = X_valid, label = y_valid)
      
      # Parameters
      params <- list(objective = "binary:logistic", eval_metric = "auc", 
      max_depth = 6, eta = 0.05, subsample = 0.8,colsample_bytree = 0.8,                  scale_pos_weight = scale_pos_weight, seed = 123)
      
      # Train with early stopping
      fit_xgb <- xgb.train(params = params, data = dtrain, nrounds = 500, 
      watchlist = list(train=dtrain, eval=dvalid),early_stopping_rounds = 50,
      verbose = 0)
      # Predict and compute AUC
      p_hat_xgb <- predict(fit_xgb, X_valid)
      roc_obj <- roc(response = valid$y_next, predictor = p_hat_xgb, quiet = TRUE)
      as.numeric(auc(roc_obj))}))

# Check performance
folds %>% summarise(mean_auc_glm = mean(auc_glm, na.rm = TRUE))
folds2 %>% summarise(mean_auc_xgb = mean(auc_xgb, na.rm = TRUE))


# filtering for more efficient modeling (however I did perform AUC with the unfiltered version as well to show the other chart)
wildfires_filtered <- wildfires_final %>% filter(FIRE_SIZE >= 300)

# proportion of positives in the dataset for comparison
baseline_pr <- mean(wildfires_filtered$y_next == 1)

evaluate_fold <- function(tr_end, va_year){
  train <- wildfires_filtered %>% filter(FIRE_YEAR <= tr_end) %>% drop_na(all_of(all.vars(frmla)))
  valid <- wildfires_filtered %>% filter(FIRE_YEAR == va_year, as.character(FIPS_CODE) %in% as.character(train$FIPS_CODE),
           STATE %in% train$STATE) %>% drop_na(all_of(all.vars(frmla)))

  # Convert response to factor
  train$y_next <- factor(train$y_next)
  valid$y_next <- factor(valid$y_next)

  # Check if train or valid has only one class; if only one class, return or it will crash
  if(length(unique(train$y_next)) < 2 || length(unique(valid$y_next)) < 2){
    return(tibble(auc_glm = NA_real_, auc_rf = NA_real_, pr_auc_glm = NA_real_, pr_auc_rf = NA_real_))}

  # Logistic Regression
  glm_fit <- glm(frmla, data = train, family = binomial)
  p_glm <- predict(glm_fit, newdata = valid, type = "response")
  df_glm <- tibble(y = valid$y_next, p = p_glm) %>% drop_na()
  auc_glm <- as.numeric(pROC::auc(pROC::roc(df_glm$y, df_glm$p, quiet=TRUE)))
  pr_auc_glm <- PRROC::pr.curve(scores.class0 = df_glm$p[df_glm$y==1], scores.class1 = df_glm$p[df_glm$y==0]
  )$auc.integral

  # Random Forest
  rf_fit <- randomForest(frmla, data = train)
  p_rf <- predict(rf_fit, newdata = valid, type="prob")[,2]
  df_rf <- tibble(y = valid$y_next, p = p_rf) %>% drop_na()
  auc_rf <- as.numeric(pROC::auc(pROC::roc(df_rf$y, df_rf$p, quiet=TRUE)))
  pr_auc_rf <- PRROC::pr.curve(scores.class0 = df_rf$p[df_rf$y==1], scores.class1 = df_rf$p[df_rf$y==0]
  )$auc.integral
  #rf_importance <- list(importance(rf_fit)) only used once because it takes a long time to process

  # Return metrics
  tibble(auc_glm = auc_glm, auc_rf = auc_rf, pr_auc_glm = pr_auc_glm, pr_auc_rf = pr_auc_rf, #rf_importance = rf_importance
  )}

# Apply to folds
folds <- folds %>%
  mutate(results = map2(train_end_year, valid_year, evaluate_fold)) %>%
  tidyr::unnest(results, names_sep = "_")

#folds$results_rf_importance[[7]]
folds %>%
  summarise(glm_auc = mean(results_auc_glm, na.rm = TRUE), rf_auc = mean(results_auc_rf,  na.rm = TRUE),
            glm_pr = mean(results_pr_auc_glm, na.rm = TRUE), rf_pr = mean(results_pr_auc_rf,  na.rm = TRUE))

# Do XGBoost separately as its more complex
features <- c("temp_3mo","ppt_3mo","FIRE_SIZE_3mo","LAT","LON", "STATE","fire_lag12", "fire_roll3", "fire_roll12",
              "size_roll3","size_roll12","MONTH")
target <- "y_next"

# Encode categorical features consistently
encode_features <- function(df, train_levels = NULL) {
  df <- df %>% mutate(across(c("STATE","MONTH"), as.factor))
  if(!is.null(train_levels)){
    for(col in c("STATE","MONTH")){
      df[[col]] <- factor(df[[col]], levels = train_levels[[col]])
      df[[col]] <- as.integer(df[[col]])}
  } else {
    for(col in c("STATE","MONTH")){
      df[[col]] <- as.integer(df[[col]])}}
  df}

# Fold setup
folds2 <- tibble(
  train_end_year = years_sorted[(min_train_years):(length(years_sorted) - 2)],
  valid_year = years_sorted[(min_train_years + 1):(length(years_sorted) - 1)]
)

folds2 <- folds2 %>%
  mutate(
    xgb_metrics = map2(train_end_year, valid_year, function(tr_end, va_year) {
      train <- wildfires_filtered %>%
        filter(FIRE_YEAR <= tr_end) %>%
        drop_na(all_of(features), all_of(target))
      
      valid <- wildfires_filtered %>%
        filter(FIRE_YEAR == va_year) %>%
        drop_na(all_of(features), all_of(target))
      
      # Skip fold if train/valid has only 1 class
      if(length(unique(train[[target]])) < 2 || length(unique(valid[[target]])) < 2){
        return(tibble(roc_auc = NA_real_, pr_auc = NA_real_))}
      
      # Save factor levels
      train_levels <- lapply(train[,c("STATE","MONTH")], levels)
      
      # Encode categorical features as integer codes
      train <- encode_features(train)
      valid <- encode_features(valid, train_levels)
      
      # Labels
      y_train <- as.numeric(train[[target]]) - 1
      y_valid <- as.numeric(valid[[target]]) - 1
      
      # Matrix for XGBoost
      X_train <- as.matrix(train[,features])
      X_valid <- as.matrix(valid[,features])
      
      # Compute scale_pos_weight for imbalance
      scale_pos_weight <- sum(y_train == 0) / sum(y_train == 1)
      
      # DMatrix
      dtrain <- xgb.DMatrix(data = X_train, label = y_train)
      dvalid <- xgb.DMatrix(data = X_valid, label = y_valid)
      
      # XGBoost parameters
      params <- list(
        objective = "binary:logistic", eval_metric = "auc", max_depth = 6, eta = 0.05,
        subsample = 0.8, colsample_bytree = 0.8, scale_pos_weight = scale_pos_weight, seed = 123)
      
      # Train XGBoost
      fit_xgb <- xgb.train(params = params, data = dtrain, nrounds = 500, watchlist = list(train = dtrain, eval = dvalid),
                           early_stopping_rounds = 50, verbose = 0)
      p_hat_xgb <- predict(fit_xgb, X_valid)
      
      # ROC-AUC & PR-AUC
      roc_obj <- pROC::roc(y_valid, p_hat_xgb, quiet = TRUE)
      roc_auc <- as.numeric(pROC::auc(roc_obj))
      pr_auc <- PRROC::pr.curve(
        scores.class0 = p_hat_xgb[y_valid == 1],
        scores.class1 = p_hat_xgb[y_valid == 0]
      )$auc.integral
      tibble(roc_auc = roc_auc, pr_auc = pr_auc)
    })
  ) %>%
  tidyr::unnest(xgb_metrics)

summary_folds2 <- folds2 %>%
  summarise(mean_roc_auc = mean(roc_auc, na.rm = TRUE), mean_pr_auc  = mean(pr_auc, na.rm = TRUE))

summary_folds1 <- folds %>%
  summarise(glm_auc = mean(results_auc_glm, na.rm = TRUE), rf_auc = mean(results_auc_rf,  na.rm = TRUE),
    glm_pr  = mean(results_pr_auc_glm, na.rm = TRUE), rf_pr = mean(results_pr_auc_rf,  na.rm = TRUE))

# Merge
combined_summary <- bind_cols(summary_folds1, summary_folds2, baseline_pr)
folds
folds2
# AUC Plot
auc_long <- folds %>%
  dplyr::select(valid_year, auc_glm, auc_rf) %>%
  pivot_longer(cols = c(auc_glm, auc_rf), names_to = "model", values_to = "AUC") %>%
  mutate(model = case_when(model == "auc_glm" ~ "GLM", model == "auc_rf"  ~ "RF")) %>%
  bind_rows(folds2 %>% dplyr::select(valid_year, auc_xgb = roc_auc) %>%
  rename(AUC = auc_xgb) %>%
  mutate(model = "XGBoost"))

AUC_plot_notable <- ggplot(auc_long, aes(x = valid_year, y = AUC, color = model)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0.4, 1)) +
  labs(title = "Model Performance (AUC) Over Validation Years (Filtered)", x = "Validation Year",
    y = "AUC", color = "Model")

