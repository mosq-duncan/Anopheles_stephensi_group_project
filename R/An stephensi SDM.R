
library(greta)
library(greta.dynamics)
library(tidyverse)
library(terra)
library(tidyterra)
library(bayesplot)
library(sf)

#load all functions
function_files <- list.files("functions/", pattern = "\\.R$", full.names = TRUE)
lapply(function_files, source)

# load hexes and hex raster lookup (all cells and populated areas)
hexes <- readRDS("output/hexes.RDS")
hex_lookup <- rast("output/rasters/derived/hex_raster_lookup.tif")
hex_populated_lookup <- rast("output/rasters/derived/hex_raster_populated_lookup.tif")

#Load covariates _ world
larval_covs <- rast("output/rasters/derived/larval_habitat_covariates.tif")
climatic_rel_abund <- rast("output/rasters/derived/climatic_rel_abund.tif")
first_detection <- readRDS("output/tabular/EA_first_detection.RDS")
mask <- rast("output/rasters/derived/mask.tif")

plot(hex_lookup)

#EA_covariates
EA_larval_covs<-resizetoEA(larval_covs)
EA_climatic_rel_abund<-resizetoEA(climatic_rel_abund)
EA_mask<-resizetoEA2(mask)

#Resample 
EA_larval_covs<- resample(EA_larval_covs, hex_lookup, method="bilinear")
EA_climatic_rel_abund<- resample(EA_climatic_rel_abund, hex_lookup, method="bilinear")
EA_mask<- resample(EA_mask, hex_lookup, method="bilinear")

# remask the hexes and everything else
EA_larval_covs <- mask(EA_larval_covs, EA_mask)
EA_climatic_rel_abund <- mask(EA_climatic_rel_abund, EA_mask)


hex_lookup <- mask(hex_lookup, EA_mask)
hex_populated_lookup <- mask(hex_populated_lookup, EA_mask)
all_hexes <- unique(hex_lookup[][, 1])
hexes <- hexes %>%
  filter(id %in% all_hexes)

# pad and a add an epsilon to the climatic relative abundance to avoid 0
# abundances later
EA_climatic_rel_abund[is.na(EA_climatic_rel_abund)] <- 0
EA_climatic_rel_abund <- mask(EA_climatic_rel_abund, EA_mask)
EA_climate_suitable <- EA_climatic_rel_abund > 0
epsilon <- .Machine$double.eps
EA_climatic_rel_abund <- EA_climatic_rel_abund + epsilon
plot(EA_climatic_rel_abund)


# filter this to points within the study region, and within the climatic
# suitability prediction
EA_climate_vals <- terra::extract(EA_climatic_rel_abund, first_detection[, 1:2])[, 2]


# augment this with background points
non_detection <- dismo::randomPoints(
  mask = raster::raster(hex_populated_lookup),
  n = 132, #changed from 500 to match detection
  p = as.matrix(first_detection[, 1:2])) %>%
  as_tibble() %>%
  mutate(
    native = FALSE,
    year_first_detected = max(first_detection$year_first_detected),
    ever_detected = 0,
    background = 1
  )

# combine these, noting which are background points (lower detection
# probability), rather than absence records
detections <- non_detection %>%
  bind_rows(
    mutate(first_detection,
           background = 0)
  ) %>%
  mutate(
    id = terra::extract(hex_populated_lookup,
                        .[, 1:2],
                        ID = FALSE)[, 1]
  )

dim(detections)

# extract the values of the larval covariates at the points
larval_covs_points <- terra::extract(EA_larval_covs,
                                     detections[, 1:2],
                                     ID = FALSE) %>% as_tibble()

plot(EA_larval_covs$built_volume)

str(larval_covs_points)



# and the climatic model
climatic_rel_abund_points <- terra::extract(EA_climatic_rel_abund,
                                            detections[, 1:2],
                                            ID = FALSE) %>% 
  as_tibble() %>%
  pull(mean)


#Combine columns
EA_An_stephensi <- detections %>% cbind(climatic_rel_abund_points, larval_covs_points)

write.csv(EA_An_stephensi, "output/EA_An_stephensi.csv")


#mantain model columns and remove NA 

EA_An_stephensi <-  EA_An_stephensi %>% 
  select(x,y,ever_detected, climatic_rel_abund_points, built_volume,tcb,tcw) %>% 
  na.omit(EA_An_stephensi)

dim(EA_An_stephensi)

#Get detection

Presence <- EA_An_stephensi$ever_detected


#Model definition of model matrix
larval_hab_formula <- ~ -1 + built_volume + tcb + tcw + climatic_rel_abund_points

Steph_model_matrix <- model.matrix(larval_hab_formula,
                                      data = EA_An_stephensi)

Steph_model_matrix<- as_data(Steph_model_matrix)
n_coefs <- ncol(Steph_model_matrix)
coefs <- normal(0, 1, dim = n_coefs)

# get log densities for points (logit)
log_larval_hab_points <- Steph_model_matrix %*% coefs

# Convert logit to probability using the logistic function
p <- ilogit(log_larval_hab_points)

# Define the likelihood: Presence/Absence follows a Bernoulli distribution
distribution(Presence) <- greta::bernoulli(p)

# Sample from the posterior using MCMC
m <- model(coefs)
plot(m)

draws <- mcmc(m, 
              n_samples = 1000, 
              warmup = 200,
              chains = 8)


bayesplot::mcmc_trace(draws)
coda::gelman.diag(draws, autoburnin = FALSE)


summary(draws)

plot(draws)

predicted_prob <- greta::calculate(p,  values = draws)

predicted_prob <- as.matrix(predicted_prob)

# Get mean probability of presence
EA_An_stephensi$predicted_presence <- colMeans(predicted_prob)

#full scale EA model
#covariates
Larva <- EA_larval_covs[[c("tcw","tcb","built_volume")]]

#Add climate abundance to the raster stack
predictor_stack <- c(Larva,EA_climatic_rel_abund)

names(predictor_stack)

Te <- predictor_stack[[4]]
names(Te)

names(predictor_stack)[4] <- "climatic_rel_abund_points"


# Convert raster stack to a data frame
predictor_df <- as.data.frame(predictor_stack, xy = TRUE, na.rm = TRUE)

# Create model matrix using the same formula as training data
larval_hab_formula <- ~ -1 + built_volume + tcb + tcw + climatic_rel_abund_points

X_new <- model.matrix(larval_hab_formula, data = predictor_df)

# Convert to greta array
X_new <- as_data(X_new)

# Compute the logit probability using trained coefficients
logit_p_new <- X_new %*% coefs
p_new <- ilogit(logit_p_new)  

# Generate predictions using greta::calculate()
predicted_prob_new <- greta::calculate(p_new, values = draws)


# Predict the posterior mean and variance from the posterior samples

predicted_prob_new_mean <- colMeans(predicted_prob_new[[1]])
predicted_prob_new_var <- apply(predicted_prob_new[[1]], 2, var)

#Insert into raster. Can use the same process for variance 
pred_rast <-  predictor_stack$tcw * 0    #make it o value raster
names(pred_rast) <- "posterior_mean"
not_na_idx <- which(!is.na(values(pred_rast)))
pred_rast[not_na_idx] <- predicted_prob_new_mean

# Save the output as a TIFF file
writeRaster(pred_rast, "output/Raster results/Anopheles_presence_probability.tif")

#Plot

ggplot() +  
  geom_spatraster(data=pred_rast) +
  geom_sf(data=EA2, color = "black", fill = NA, linewidth = 0.5) + 
  scale_fill_viridis_c() +  
  theme_minimal()






