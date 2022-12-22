# Wolbachia TPP economic models
# combines epidemiological and economic burden data to identify which areas are 
# most economically efficient to treat reduce the global burden of dengue by
# fixed percentages and the cost thresholds identified by these thresholds

# Author: Oliver Brady (oliver.brady@lshtm.ac.uk)
# Date: 08 Dec 2022


rm(list = ls())

# packages
require(maps)
require(raster)
require(sf)
require(tmap)
require(exactextractr)
require(ggplot2)
require(ggrepel)

# set working directory
setwd("/Users/eideobra/Dropbox/07_Indonesia/14_TPP")


## 01 load in datasets ###

# load in dengue burden map (Bhatt et al. Nature 2013) and set relevant CRS
denBurden <- raster("TPP_data/Dengue_burden_symp1508.asc")
denBurden_low <- raster("TPP_data/Dengue_burden_symp1508_low.asc")
denBurden_high <- raster("TPP_data/Dengue_burden_symp1508_high.asc")
crs(denBurden) <- crs(denBurden_low) <- crs(denBurden_high) <- 
  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# population data
pop_5km <- raster("TPP_data/Population_raster_5km.tiff")
crs(pop_5km) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# load in country shapefiles
ad0 <- read_sf("TPP_data/Admin0(2011)/admin0.shp")
ad2 <- read_sf("TPP_data/Admin2(2011)/admin2.shp")
ad0_ras <- raster("TPP_data/Admin0(2011)/admin0_5k_raster.tif") # admin0 in raster format

# load in average cost per symptomatic case raster from the Shepard data
cost <- raster("TPP_data/Shepard_Avg_Dir_med_cost_per_symp_case.tif")

# load in cost of routine vector control (per capita per year) raster
routineVC_cost = raster("TPP_data/Routine_vector_control_costs_per_capita_per_year.tif")

# gdp per capita data
GDP_percap <- read.csv("TPP_data/WB_GDP_2020_data.csv")

# GDP deflators
GDP_deflator_tab <- read.csv("TPP_data/WB_GDP_deflator.csv")


# function to adjust to 2020 USD using the GDPdeflators
GDP.deflator <- function(dataYear, targetYear){
  
  # match target columns
  colDY = match(dataYear, gsub("X", "", colnames(GDP_deflator_tab)))
  colTY = match(targetYear, gsub("X", "", colnames(GDP_deflator_tab)))
  
  # assemble records
  DY_records =GDP_deflator_tab[, colDY]
  TY_records =GDP_deflator_tab[, colTY]
  
  # calculate national multiplier to go from data year to target year
  Mult <- TY_records / DY_records
  
  # fill in NAs (small number of countries missin values) with the median value
  Mult[is.na(Mult)] = median(Mult, na.rm = T)
  
  return(data.frame(ISO3 = GDP_deflator_tab$Country.Code,
                    Multiplier = Mult))
}




#### 02 Calculating effectiveness and averted cases ###

Effectiveness <- 0.7 # universal target from the TPP

# post Wolbachia dengue burden and averted cases
PW_denBurden_low <- denBurden_low * (1 - Effectiveness)
PW_denBurden_mid <- denBurden * (1 - Effectiveness)
PW_denBurden_high <- denBurden_high * (1 - Effectiveness)

PW_AvertedCases_low <- denBurden_low - PW_denBurden_low
PW_AvertedCases_mid <- denBurden - PW_denBurden_mid
PW_AvertedCases_high <- denBurden_high - PW_denBurden_high




### 03 Direct medical benefits and avertable vector control costs ###

# percentage of routine vector control costs that are avertable (i.e. the emergency costs) -only 3 month duration
Averted_VC_cost_low = routineVC_cost * pop_5km * 0.20 * (3 / 12)
Averted_VC_cost_mid = routineVC_cost * pop_5km * 0.35 * (3 / 12)
Averted_VC_cost_high = routineVC_cost * pop_5km * 0.50 * (3 / 12)

# avertable insecticide costs (Optimal TPP category)
Averted_Insecticide_cost_low = Averted_VC_cost_low * 0.01
Averted_Insecticide_cost_mid = Averted_VC_cost_mid * 0.10
Averted_Insecticide_cost_high = Averted_VC_cost_high * 0.20

PW_AnnualSavings_low <- PW_AvertedCases_low * cost + Averted_VC_cost_low
PW_AnnualSavings_mid <- PW_AvertedCases_mid * cost + Averted_VC_cost_mid
PW_AnnualSavings_high <- PW_AvertedCases_high * cost + Averted_VC_cost_high

PW_AnnualSavings_Preferred_low <- Averted_VC_cost_low
PW_AnnualSavings_Preferred_mid <- Averted_VC_cost_mid
PW_AnnualSavings_Preferred_high <- Averted_VC_cost_high

PW_AnnualSavings_Optimal_low <- Averted_Insecticide_cost_low
PW_AnnualSavings_Optimal_mid <- Averted_Insecticide_cost_mid
PW_AnnualSavings_Optimal_high <- Averted_Insecticide_cost_high

# work out cost discounted multipliers for different time horizons (3% per year)
M_10 <- sum(0.97^(0:(10 - 1)))
M_5 <- sum(0.97^(0:(5 - 1)))
M_3 <- sum(0.97^(0:(3 - 1)))

PW_Savings_low_10yr = PW_AnnualSavings_low * M_10
PW_Savings_low_5yr = PW_AnnualSavings_low * M_5
PW_Savings_low_3yr = PW_AnnualSavings_low * M_3

PW_Savings_mid_10yr = PW_AnnualSavings_mid * M_10
PW_Savings_mid_5yr = PW_AnnualSavings_mid * M_5
PW_Savings_mid_3yr = PW_AnnualSavings_mid * M_3

PW_Savings_high_10yr = PW_AnnualSavings_high * M_10
PW_Savings_high_5yr = PW_AnnualSavings_high * M_5
PW_Savings_high_3yr = PW_AnnualSavings_high * M_3


PW_Savings_Preferred_low_10yr = PW_AnnualSavings_Preferred_low * M_10
PW_Savings_Preferred_low_5yr = PW_AnnualSavings_Preferred_low * M_5
PW_Savings_Preferred_low_3yr = PW_AnnualSavings_Preferred_low * M_3

PW_Savings_Preferred_mid_10yr = PW_AnnualSavings_Preferred_mid * M_10
PW_Savings_Preferred_mid_5yr = PW_AnnualSavings_Preferred_mid * M_5
PW_Savings_Preferred_mid_3yr = PW_AnnualSavings_Preferred_mid * M_3

PW_Savings_Preferred_high_10yr = PW_AnnualSavings_Preferred_high * M_10
PW_Savings_Preferred_high_5yr = PW_AnnualSavings_Preferred_high * M_5
PW_Savings_Preferred_high_3yr = PW_AnnualSavings_Preferred_high * M_3


PW_Savings_Optimal_low_10yr = PW_AnnualSavings_Optimal_low * M_10
PW_Savings_Optimal_low_5yr = PW_AnnualSavings_Optimal_low * M_5
PW_Savings_Optimal_low_3yr = PW_AnnualSavings_Optimal_low * M_3

PW_Savings_Optimal_mid_10yr = PW_AnnualSavings_Optimal_mid * M_10
PW_Savings_Optimal_mid_5yr = PW_AnnualSavings_Optimal_mid * M_5
PW_Savings_Optimal_mid_3yr = PW_AnnualSavings_Optimal_mid * M_3

PW_Savings_Optimal_high_10yr = PW_AnnualSavings_Optimal_high * M_10
PW_Savings_Optimal_high_5yr = PW_AnnualSavings_Optimal_high * M_5
PW_Savings_Optimal_high_3yr = PW_AnnualSavings_Optimal_high * M_3




### 04 costs of Wolbachia ###

# load wolbachia costing models from Brady et al. 2020
load("TPP_data/Cost_model.RData")

# vectorise population raster to make it easier to work with
pop_5km_v <- as.vector(pop_5km)

# corresponding vector of GDPs
ad0_ras_v = as.vector(ad0_ras)
ad0_ras_v = ad0$COUNTRY_ID[match(ad0_ras_v, ad0$GAUL_CODE)]
GDP_ras_v = GDP_percap$GDP_percap_2020[match(ad0_ras_v, GDP_percap$Country.Code)]
#for areas where we do have population data, but not GDP data, assume median GDP per capita
GDP_ras_v[(is.na(GDP_ras_v) & (!is.na(pop_5km_v)))] = median(GDP_ras_v, na.rm = T)


# predict Phase 1 and Phase 2 Wolbachia programme costs from the Brady et al. cost model
P1P2C <- predict(glmmod,
                 newdata = data.frame(Pop_density =  log(pop_5km_v / 25, 10), # as 5k x 5km cells
                                      Programme.phase = 1,
                                      Type_release = "EGGS",
                                      GDP = GDP_ras_v),
                 se.fit = F,
                 type = "response")

P2C <- predict(glmmod,
               newdata = data.frame(Pop_density =  log(pop_5km_v / 25, 10),
                                    Programme.phase = 2,
                                    Type_release = "EGGS",
                                    GDP = GDP_ras_v),
               se.fit = F,
               type = "response")

# transform back to response scale and scale up to 5 x 5km grid
P1P2C = (10 ^ P1P2C) * 5^2
P2C = (10 ^ P2C) * 5^2

# phase 3 and 4 (monitoring) are an pre-agreed proportion of phase 2 costs (see Brady et al. paper)
P1P2P3P4C <- P1P2C + P2C * 0.16 + P2C * 0.08

# adjust between 2018 and 2020 USD
Multiplier_tab <- GDP.deflator(2018, 2020)
Multiplier_ras_v <- Multiplier_tab$Multiplier[match(ad0_ras_v, Multiplier_tab$ISO3)]
# assign median value for areas without GDP data
Multiplier_ras_v[is.na(Multiplier_ras_v)] = median(Multiplier_ras_v, na.rm = T)

# calcualte combined costs of all phases of the Wolbachia programme
P1P2C = P1P2C * Multiplier_ras_v
P2C = P2C * Multiplier_ras_v
P1P2P3P4C = P1P2P3P4C * Multiplier_ras_v


# total costs of different business models put back into rasters
TurnKey_costs <- pop_5km
IS_costs_mid <- IS_costs_low <- IS_costs_high <- pop_5km
values(TurnKey_costs) <- P1P2P3P4C

# proportion from budget dissection
Production_proportions <- c(18.6267219,	10.12876394,	44.29571442)
values(IS_costs_mid) <- (P1P2C) * Production_proportions[1]
values(IS_costs_low) <- (P1P2C) * Production_proportions[2]
values(IS_costs_high) <- (P1P2C) * Production_proportions[3]







### 05 pixel level analyses ####

#compile 5 x 5km pixels into a big dataframe
world.pixels <- data.frame(ID = 1:length(pop_5km_v),
                           GAUL = as.vector(ad0_ras),
                           Pop = pop_5km_v,
                           Averted_cases_low = as.vector(PW_AvertedCases_low),
                           Averted_cases_mid = as.vector(PW_AvertedCases_mid),
                           Averted_cases_high = as.vector(PW_AvertedCases_high),
                           Savings_10yr_low = as.vector(PW_Savings_low_10yr),
                           Savings_10yr_mid = as.vector(PW_Savings_mid_10yr),
                           Savings_10yr_high = as.vector(PW_Savings_high_10yr),
                           Savings_5yr_low = as.vector(PW_Savings_low_5yr),
                           Savings_5yr_mid = as.vector(PW_Savings_mid_5yr),
                           Savings_5yr_high = as.vector(PW_Savings_high_5yr),
                           Savings_3yr_low = as.vector(PW_Savings_low_3yr),
                           Savings_3yr_mid = as.vector(PW_Savings_mid_3yr),
                           Savings_3yr_high = as.vector(PW_Savings_high_3yr),
                           Savings_Preferred_10yr_low = as.vector(PW_Savings_Preferred_low_10yr),
                           Savings_Preferred_10yr_mid = as.vector(PW_Savings_Preferred_mid_10yr),
                           Savings_Preferred_10yr_high = as.vector(PW_Savings_Preferred_high_10yr),
                           Savings_Preferred_5yr_low = as.vector(PW_Savings_Preferred_low_5yr),
                           Savings_Preferred_5yr_mid = as.vector(PW_Savings_Preferred_mid_5yr),
                           Savings_Preferred_5yr_high = as.vector(PW_Savings_Preferred_high_5yr),
                           Savings_Preferred_3yr_low = as.vector(PW_Savings_Preferred_low_3yr),
                           Savings_Preferred_3yr_mid = as.vector(PW_Savings_Preferred_mid_3yr),
                           Savings_Preferred_3yr_high = as.vector(PW_Savings_Preferred_high_3yr),
                           Savings_Optimal_10yr_low = as.vector(PW_Savings_Optimal_low_10yr),
                           Savings_Optimal_10yr_mid = as.vector(PW_Savings_Optimal_mid_10yr),
                           Savings_Optimal_10yr_high = as.vector(PW_Savings_Optimal_high_10yr),
                           Savings_Optimal_5yr_low = as.vector(PW_Savings_Optimal_low_5yr),
                           Savings_Optimal_5yr_mid = as.vector(PW_Savings_Optimal_mid_5yr),
                           Savings_Optimal_5yr_high = as.vector(PW_Savings_Optimal_high_5yr),
                           Savings_Optimal_3yr_low = as.vector(PW_Savings_Optimal_low_3yr),
                           Savings_Optimal_3yr_mid = as.vector(PW_Savings_Optimal_mid_3yr),
                           Savings_Optimal_3yr_high = as.vector(PW_Savings_Optimal_high_3yr),
                           Turnkey_costs = P1P2P3P4C,
                           Integrated_supply_costs_low = (P1P2C) * Production_proportions[2],
                           Integrated_supply_costs_mid = (P1P2C) * Production_proportions[1],
                           Integrated_supply_costs_high = (P1P2C) * Production_proportions[3])

# trim to dengue endemic areas with > 0 population and remove those that we don't have cost data for
world.pixels = world.pixels[!is.na(world.pixels$Averted_cases_mid), ]
world.pixels = world.pixels[(world.pixels$Averted_cases_mid > 0), ]
world.pixels = world.pixels[(world.pixels$Pop > 0), ]
world.pixels = world.pixels[!is.na(world.pixels$Savings_10yr_low), ]

# calculate benefit / cost ratio for current (estimated) programme costs
world.pixels$BCR_TK <- world.pixels$Savings_10yr_mid / world.pixels$Turnkey_costs
world.pixels$BCR_TK_VC <- world.pixels$Savings_Preferred_10yr_mid / world.pixels$Turnkey_costs
world.pixels$BCR_IS_mid <- world.pixels$Savings_10yr_mid / world.pixels$Integrated_supply_costs_mid


### 06 TPP tables

##  minimum criteria:

# rank from most cost efficient to least cost efficient 
# efficient = benefit to cost ratio
world.pixels = world.pixels[order(world.pixels$BCR_TK, decreasing = T), ]


# calculate pixels that need to be included to meet WHO targets
cum_averted_cases_low = cumsum(world.pixels$Averted_cases_low)
cum_averted_cases_mid = cumsum(world.pixels$Averted_cases_mid)
cum_averted_cases_high = cumsum(world.pixels$Averted_cases_high)

gburd = sum(as.vector(denBurden), na.rm = T) # global burden

pix_cut_low = which.max(cum_averted_cases_low >= (0.25 * gburd))
pix_cut_mid = which.max(cum_averted_cases_mid >= (0.25 * gburd))
pix_cut_high = which.max(cum_averted_cases_high >= (0.25 * gburd))

#if only half the reduction in global burden needs to be achieved using Wolbachia
pix_cut50_low = which.max(cum_averted_cases_low >= (0.25 * 0.5 * gburd))
pix_cut50_mid = which.max(cum_averted_cases_mid >= (0.25 * 0.5 * gburd))
pix_cut50_high = which.max(cum_averted_cases_high >= (0.25 * 0.5 * gburd))


## for Preferred and Optimal TPP
# rank from most cost efficient to least cost efficient from Vector control avertable costs perspective
# efficient = benefit to cost ratio (just for use in the optimal category)
world.pixels_VC = world.pixels[order(world.pixels$BCR_TK_VC, decreasing = T), ]

# pixels needed to meet WHO targets
cum_averted_cases_low = cumsum(world.pixels_VC$Averted_cases_low)
cum_averted_cases_mid = cumsum(world.pixels_VC$Averted_cases_mid)
cum_averted_cases_high = cumsum(world.pixels_VC$Averted_cases_high)

pix_cut_VC_low = which.max(cum_averted_cases_low >= (0.25 * gburd))
pix_cut_VC_mid = which.max(cum_averted_cases_mid >= (0.25 * gburd))
pix_cut_VC_high = which.max(cum_averted_cases_high >= (0.25 * gburd))

#if only half the reduction in global burden needs to be achieved using Wolbachia
pix_cut50_VC_low = which.max(cum_averted_cases_low >= (0.25 * 0.5 * gburd))
pix_cut50_VC_mid = which.max(cum_averted_cases_mid >= (0.25 * 0.5 * gburd))
pix_cut50_VC_high = which.max(cum_averted_cases_high >= (0.25 * 0.5 * gburd))


# calculate long term savings (i.e. programme break even cost) / people protected
# to give target cost per person
world.pixels$TPP_10yr_low = world.pixels$Savings_10yr_low / world.pixels$Pop
world.pixels$TPP_10yr_mid = world.pixels$Savings_10yr_mid / world.pixels$Pop
world.pixels$TPP_10yr_high = world.pixels$Savings_10yr_high / world.pixels$Pop

world.pixels$TPP_5yr_low = world.pixels$Savings_5yr_low / world.pixels$Pop
world.pixels$TPP_5yr_mid = world.pixels$Savings_5yr_mid / world.pixels$Pop
world.pixels$TPP_5yr_high = world.pixels$Savings_5yr_high / world.pixels$Pop

world.pixels$TPP_3yr_low = world.pixels$Savings_3yr_low / world.pixels$Pop
world.pixels$TPP_3yr_mid = world.pixels$Savings_3yr_mid / world.pixels$Pop
world.pixels$TPP_3yr_high = world.pixels$Savings_3yr_high / world.pixels$Pop


world.pixels_VC$TPP_Preferred_10yr_low = world.pixels_VC$Savings_Preferred_10yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_10yr_mid = world.pixels_VC$Savings_Preferred_10yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_10yr_high = world.pixels_VC$Savings_Preferred_10yr_high / world.pixels_VC$Pop

world.pixels_VC$TPP_Preferred_5yr_low = world.pixels_VC$Savings_Preferred_5yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_5yr_mid = world.pixels_VC$Savings_Preferred_5yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_5yr_high = world.pixels_VC$Savings_Preferred_5yr_high / world.pixels_VC$Pop

world.pixels_VC$TPP_Preferred_3yr_low = world.pixels_VC$Savings_Preferred_3yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_3yr_mid = world.pixels_VC$Savings_Preferred_3yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Preferred_3yr_high = world.pixels_VC$Savings_Preferred_3yr_high / world.pixels_VC$Pop


world.pixels_VC$TPP_Optimal_10yr_low = world.pixels_VC$Savings_Optimal_10yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_10yr_mid = world.pixels_VC$Savings_Optimal_10yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_10yr_high = world.pixels_VC$Savings_Optimal_10yr_high / world.pixels_VC$Pop

world.pixels_VC$TPP_Optimal_5yr_low = world.pixels_VC$Savings_Optimal_5yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_5yr_mid = world.pixels_VC$Savings_Optimal_5yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_5yr_high = world.pixels_VC$Savings_Optimal_5yr_high / world.pixels_VC$Pop

world.pixels_VC$TPP_Optimal_3yr_low = world.pixels_VC$Savings_Optimal_3yr_low / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_3yr_mid = world.pixels_VC$Savings_Optimal_3yr_mid / world.pixels_VC$Pop
world.pixels_VC$TPP_Optimal_3yr_high = world.pixels_VC$Savings_Optimal_3yr_high / world.pixels_VC$Pop



# assembling TPP tables

# Minimum
# cost to reach all averas
TPP_table <- data.frame(TPP_10y_l = min(world.pixels$TPP_10yr_low[1:pix_cut_low]),
                        TPP_10y_m = min(world.pixels$TPP_10yr_mid[1:pix_cut_mid]),
                        TPP_10y_h = min(world.pixels$TPP_10yr_high[1:pix_cut_high]),
                        TPP_5y_l = min(world.pixels$TPP_5yr_low[1:pix_cut_low]),
                        TPP_5y_m = min(world.pixels$TPP_5yr_mid[1:pix_cut_mid]),
                        TPP_5y_h = min(world.pixels$TPP_5yr_high[1:pix_cut_high]),
                        TPP_3y_l = min(world.pixels$TPP_3yr_low[1:pix_cut_low]),
                        TPP_3y_m = min(world.pixels$TPP_3yr_mid[1:pix_cut_mid]),
                        TPP_3y_h = min(world.pixels$TPP_3yr_high[1:pix_cut_high]))
TPP_table <- rbind(TPP_table, data.frame(TPP_10y_l = min(world.pixels$TPP_10yr_low[1:pix_cut50_low]),
                                         TPP_10y_m = min(world.pixels$TPP_10yr_mid[1:pix_cut50_mid]),
                                         TPP_10y_h = min(world.pixels$TPP_10yr_high[1:pix_cut50_high]),
                                         TPP_5y_l = min(world.pixels$TPP_5yr_low[1:pix_cut50_low]),
                                         TPP_5y_m = min(world.pixels$TPP_5yr_mid[1:pix_cut50_mid]),
                                         TPP_5y_h = min(world.pixels$TPP_5yr_high[1:pix_cut50_high]),
                                         TPP_3y_l = min(world.pixels$TPP_3yr_low[1:pix_cut50_low]),
                                         TPP_3y_m = min(world.pixels$TPP_3yr_mid[1:pix_cut50_mid]),
                                         TPP_3y_h = min(world.pixels$TPP_3yr_high[1:pix_cut50_high])))
row.names(TPP_table) <- c("Full contribution", "50% contribution")

# view and save
TPP_table
write.csv(TPP_table, file = "TPP_outputs/TPP_minimum_summary_table.csv")

# hold pix_cut_mid for later analysis
cost_burden_peg <- pix_cut_mid

# average cost within areas (not used)
#TPP_table_avg <- data.frame(TPP_10y_l = median(world.pixels$TPP_10yr_low[1:pix_cut_low]),
#                            TPP_10y_m = median(world.pixels$TPP_10yr_mid[1:pix_cut_mid]),
#                            TPP_10y_h = median(world.pixels$TPP_10yr_high[1:pix_cut_high]),
#                            TPP_5y_l = median(world.pixels$TPP_5yr_low[1:pix_cut_low]),
#                            TPP_5y_m = median(world.pixels$TPP_5yr_mid[1:pix_cut_mid]),
#                            TPP_5y_h = median(world.pixels$TPP_5yr_high[1:pix_cut_high]),
#                            TPP_3y_l = median(world.pixels$TPP_3yr_low[1:pix_cut_low]),
#                            TPP_3y_m = median(world.pixels$TPP_3yr_mid[1:pix_cut_mid]),
#                            TPP_3y_h = median(world.pixels$TPP_3yr_high[1:pix_cut_high]))
#TPP_table_avg <- rbind(TPP_table_avg, data.frame(TPP_10y_l = median(world.pixels$TPP_10yr_low[1:pix_cut50_low]),
#                                                 TPP_10y_m = median(world.pixels$TPP_10yr_mid[1:pix_cut50_mid]),
#                                                 TPP_10y_h = median(world.pixels$TPP_10yr_high[1:pix_cut50_high]),
#                                                 TPP_5y_l = median(world.pixels$TPP_5yr_low[1:pix_cut50_low]),
#                                                 TPP_5y_m = median(world.pixels$TPP_5yr_mid[1:pix_cut50_mid]),
#                                                 TPP_5y_h = median(world.pixels$TPP_5yr_high[1:pix_cut50_high]),
#                                                 TPP_3y_l = median(world.pixels$TPP_3yr_low[1:pix_cut50_low]),
#                                                 TPP_3y_m = median(world.pixels$TPP_3yr_mid[1:pix_cut50_mid]),
#                                                 TPP_3y_h = median(world.pixels$TPP_3yr_high[1:pix_cut50_high])))
#row.names(TPP_table_avg) <- c("Full contribution", "50% contribution")
#TPP_table_avg


# overall costs (billions USD)
cost_tab = rbind(TPP_table[1, c(2, 5, 8)] * sum(world.pixels$Pop[1:pix_cut_mid]),
                 TPP_table[2, c(2, 5, 8)] * sum(world.pixels$Pop[1:pix_cut50_mid]))
round(cost_tab / 1000000000, 2)


# Preferred
TPP_table <- data.frame(TPP_10y_l = min(world.pixels_VC$TPP_Preferred_10yr_low[1:pix_cut_VC_low]),
                        TPP_10y_m = min(world.pixels_VC$TPP_Preferred_10yr_mid[1:pix_cut_VC_mid]),
                        TPP_10y_h = min(world.pixels_VC$TPP_Preferred_10yr_high[1:pix_cut_VC_high]),
                        TPP_5y_l = min(world.pixels_VC$TPP_Preferred_5yr_low[1:pix_cut_VC_low]),
                        TPP_5y_m = min(world.pixels_VC$TPP_Preferred_5yr_mid[1:pix_cut_VC_mid]),
                        TPP_5y_h = min(world.pixels_VC$TPP_Preferred_5yr_high[1:pix_cut_VC_high]),
                        TPP_3y_l = min(world.pixels_VC$TPP_Preferred_3yr_low[1:pix_cut_VC_low]),
                        TPP_3y_m = min(world.pixels_VC$TPP_Preferred_3yr_mid[1:pix_cut_VC_mid]),
                        TPP_3y_h = min(world.pixels_VC$TPP_Preferred_3yr_high[1:pix_cut_VC_high]))
TPP_table <- rbind(TPP_table, data.frame(TPP_10y_l = min(world.pixels_VC$TPP_Preferred_10yr_low[1:pix_cut50_VC_low]),
                                         TPP_10y_m = min(world.pixels_VC$TPP_Preferred_10yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_10y_h = min(world.pixels_VC$TPP_Preferred_10yr_high[1:pix_cut50_VC_high]),
                                         TPP_5y_l = min(world.pixels_VC$TPP_Preferred_5yr_low[1:pix_cut50_VC_low]),
                                         TPP_5y_m = min(world.pixels_VC$TPP_Preferred_5yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_5y_h = min(world.pixels_VC$TPP_Preferred_5yr_high[1:pix_cut50_VC_high]),
                                         TPP_3y_l = min(world.pixels_VC$TPP_Preferred_3yr_low[1:pix_cut50_VC_low]),
                                         TPP_3y_m = min(world.pixels_VC$TPP_Preferred_3yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_3y_h = min(world.pixels_VC$TPP_Preferred_3yr_high[1:pix_cut50_VC_high])))
row.names(TPP_table) <- c("Full contribution", "50% contribution")

# view and save
TPP_table
write.csv(TPP_table, file = "TPP_outputs/TPP_preferred_summary_table.csv")

# small modification for production costs estimate in TPP - preferred table
(Production_proportions[c(2,1,3)] / 100) * TPP_table[1, 7:9]



# Optimal
TPP_table <- data.frame(TPP_10y_l = min(world.pixels_VC$TPP_Optimal_10yr_low[1:pix_cut_VC_low]),
                        TPP_10y_m = min(world.pixels_VC$TPP_Optimal_10yr_mid[1:pix_cut_VC_mid]),
                        TPP_10y_h = min(world.pixels_VC$TPP_Optimal_10yr_high[1:pix_cut_VC_high]),
                        TPP_5y_l = min(world.pixels_VC$TPP_Optimal_5yr_low[1:pix_cut_VC_low]),
                        TPP_5y_m = min(world.pixels_VC$TPP_Optimal_5yr_mid[1:pix_cut_VC_mid]),
                        TPP_5y_h = min(world.pixels_VC$TPP_Optimal_5yr_high[1:pix_cut_VC_high]),
                        TPP_3y_l = min(world.pixels_VC$TPP_Optimal_3yr_low[1:pix_cut_VC_low]),
                        TPP_3y_m = min(world.pixels_VC$TPP_Optimal_3yr_mid[1:pix_cut_VC_mid]),
                        TPP_3y_h = min(world.pixels_VC$TPP_Optimal_3yr_high[1:pix_cut_VC_high]))
TPP_table <- rbind(TPP_table, data.frame(TPP_10y_l = min(world.pixels_VC$TPP_Optimal_10yr_low[1:pix_cut50_VC_low]),
                                         TPP_10y_m = min(world.pixels_VC$TPP_Optimal_10yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_10y_h = min(world.pixels_VC$TPP_Optimal_10yr_high[1:pix_cut50_VC_high]),
                                         TPP_5y_l = min(world.pixels_VC$TPP_Optimal_5yr_low[1:pix_cut50_VC_low]),
                                         TPP_5y_m = min(world.pixels_VC$TPP_Optimal_5yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_5y_h = min(world.pixels_VC$TPP_Optimal_5yr_high[1:pix_cut50_VC_high]),
                                         TPP_3y_l = min(world.pixels_VC$TPP_Optimal_3yr_low[1:pix_cut50_VC_low]),
                                         TPP_3y_m = min(world.pixels_VC$TPP_Optimal_3yr_mid[1:pix_cut50_VC_mid]),
                                         TPP_3y_h = min(world.pixels_VC$TPP_Optimal_3yr_high[1:pix_cut50_VC_high])))
row.names(TPP_table) <- c("Full contribution", "50% contribution")

# view and save
TPP_table
write.csv(TPP_table, file = "TPP_outputs/TPP_optimal_summary_table.csv")









### 07 maps of which areas would be treated ###

# assign treated yes/no to rasters
baselayer = ad0_ras >= 0
ID_vec = 1:length(pop_5km_v)
ID_vec = ID_vec %in% world.pixels$ID[1:pix_cut_mid]
ID_vec[is.na(as.vector(baselayer))] = NA
examras1 = examras1b = baselayer
values(examras1) = (values(examras1) + ID_vec)

# and map of their costs (for extraction purposes only)
global_cost_vec_100 = (world.pixels$Savings_5yr_mid / world.pixels$Pop)[match(1:length(pop_5km_v), world.pixels$ID[1:pix_cut_mid])]
global_cost_vec_100[is.na(global_cost_vec_100)] = 9999999
values(examras1b) = global_cost_vec_100

# simple map of targetted areas
p1 <- tm_shape(examras1) +
  tm_raster(legend.show = FALSE)

# save
tmap_save(p1, file = "TPP_outputs/Treat_areas_100perc.tiff")
writeRaster(examras1, file = "TPP_outputs/Treat_areas_100perc_raster.tiff")


# repeat process for areas targetted for a 12.5% global burden reduction
ID_vec50 = 1:length(pop_5km_v)
ID_vec50 = ID_vec50 %in% world.pixels$ID[1:pix_cut50_mid]
ID_vec50[is.na(as.vector(baselayer))] = NA
examras50 = examras50b = baselayer
values(examras50) = (values(examras50) + ID_vec50)

# and map of their costs (for extraction purposes only)
global_cost_vec_50 = (world.pixels$Savings_5yr_mid / world.pixels$Pop)[match(1:length(pop_5km_v), world.pixels$ID[1:pix_cut50_mid])]
global_cost_vec_50[is.na(global_cost_vec_50)] = 9999999
values(examras50b) = global_cost_vec_50

# simple map of targetted areas
p2 <- tm_shape(examras50) +
  tm_raster(legend.show = FALSE)

# save
tmap_save(p2, file = "TPP_outputs/Treat_areas_50perc.tiff")
writeRaster(examras50, file = "TPP_outputs/Treat_areas_50perc_raster.tiff")




### 08 municipality-level summaries ###

# convert raster to logical treated / not treated
examras1  = examras1 == 2

# area extraction
sf::sf_use_s2(FALSE)
mun_area <- as.numeric(st_area(ad2)) / (1000^2)
prop_covered <- exact_extract(examras1, ad2, fun = "sum") / exact_extract(examras1^0, ad2, fun = "sum")

# draw together into one big dataframe
mun_sum <- data.frame(A2_NAME = ad2$NAME, # municipality name
                      A2_GAUL_CODE = ad2$GAUL_CODE, # municipality GAUL code
                      A2COUNTRY = ad2$COUNTRY_ID, # three letter ISO country name
                      Area_mun = mun_area, # area of the municipality (Km2)
                      Area_covered = mun_area * prop_covered, # area of the municipality to be covered with Wolbachia (Km2)
                      Pop_mun = exact_extract(pop_5km, ad2, fun = "sum"), # total population of the municipality
                      Pop_covered = exact_extract(pop_5km * examras1, ad2, fun = "sum"), # total population of the municipality covered
                      Burden_mun = exact_extract(denBurden, ad2, fun = "sum"), # estimated dengue burden of the municipality
                      Burden_covered = exact_extract(denBurden * examras1, ad2, fun = "sum"), # burden of the municipality in areas where Wolbachia will be released
                      Target_cost = exact_extract(examras1b, ad2, fun = "min"), # cost at which Wolbachia can be implemented and meet TPP minimum criteria
                      VC_cost_per_year_covered = exact_extract(routineVC_cost * pop_5km * examras1, ad2, fun = "sum")) # estimated total vector control spend in the municipality 
# remove untreated areas
mun_sum = mun_sum[mun_sum$Area_covered > 0, ]
# save file
write.csv(mun_sum, file = "TPP_outputs/Munsum_global.csv")



# municipality-level summaries - 12.5% global burden reduction
examras50  = examras50 == 2

# area extraction
prop_covered50 <- exact_extract(examras50, ad2, fun = "sum") / exact_extract(examras50^0, ad2, fun = "sum")

# variable coding the same as before
mun_sum50 <- data.frame(A2_NAME = ad2$NAME,
                        A2_GAUL_CODE = ad2$GAUL_CODE,
                        A2COUNTRY = ad2$COUNTRY_ID,
                        Area_mun = mun_area,
                        Area_covered = mun_area * prop_covered50,
                        Pop_mun = exact_extract(pop_5km, ad2, fun = "sum"),
                        Pop_covered = exact_extract(pop_5km * examras50, ad2, fun = "sum"),
                        Burden_mun = exact_extract(denBurden, ad2, fun = "sum"),
                        Burden_covered = exact_extract(denBurden * examras50, ad2, fun = "sum"),
                        Target_cost = exact_extract(examras50b, ad2, fun = "min"),
                        VC_cost_per_year_covered = exact_extract(routineVC_cost * pop_5km * examras50, ad2, fun = "sum"))
# remove untreated areas
mun_sum50 = mun_sum50[mun_sum50$Area_covered > 0, ]
# save file
write.csv(mun_sum50, file = "TPP_outputs/Munsum_global_12_5perc.csv")




# quick summary of areas
# Area
area_ras = area(ad0_ras)
sum(as.vector(area_ras)[world.pixels$ID[1:pix_cut_mid]], na.rm = T)
sum(as.vector(area_ras)[world.pixels$ID[1:pix_cut50_mid]], na.rm = T)

# countries
length(unique(world.pixels$GAUL[1:pix_cut_mid]))
length(unique(world.pixels$GAUL[1:pix_cut50_mid]))

# how many areas would already meet this target under current cost base?
wp_sub <- world.pixels[1:pix_cut_mid, ]
100 * sum(wp_sub$Savings_10yr_mid >= wp_sub$Turnkey_costs) / nrow(wp_sub)
# if mosquito control costs were 50% of medical costs
100 * sum((wp_sub$Savings_10yr_mid * 1.5) >= wp_sub$Turnkey_costs) / nrow(wp_sub)

# what % of dengue at risk area?
100 * pix_cut_mid / sum(as.vector(denBurden) > 0, na.rm = T)
100 * pix_cut50_mid / sum(as.vector(denBurden) > 0, na.rm = T)

# what % of dengue at risk urban area?
treated_areas = pop_5km_v[world.pixels$ID[1:pix_cut_mid]]
quantile(treated_areas, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)) / 25
pop_areas = ((as.vector(denBurden) > 0) + (pop_5km_v > (300 * 5^2))) > 1
100 * pix_cut_mid / sum(pop_areas, na.rm = T)

treated_areas = pop_5km_v[world.pixels$ID[1:pix_cut50_mid]]
quantile(treated_areas, probs = c(0.01, 0.05, 0.5, 0.95, 0.99)) / 25
pop_areas = ((as.vector(denBurden) > 0) + (pop_5km_v > (300 * 5^2))) > 1
100 * pix_cut_mid / sum(pop_areas, na.rm = T)






#### 09 pixels needed to meet national burden reduction targets ###

# work out burden in each county
c_burd <- data.frame(GAUL = as.vector(ad0_ras),
                     Burd = as.vector(denBurden))

c_burd <- aggregate(Burd ~ GAUL, data = c_burd, FUN = sum, na.rm = T)

# collectors
u_GAULS <- unique(world.pixels$GAUL)

# trim to higher dengue burden countries (total annual symptomatic cases > 10k)?
trim =TRUE
if(trim){
  u_GAULS = u_GAULS[u_GAULS%in% c_burd$GAUL[c_burd$Burd > 10000]]
}

fc_TPP <- data.frame(GAUL = u_GAULS,
                     Y10_l =NA,
                     Y10_m =NA,
                     Y10_h =NA,
                     Y5_l =NA,
                     Y5_m =NA,
                     Y5_h =NA,
                     Y3_l =NA,
                     Y3_m =NA,
                     Y3_h =NA)
hc_TPP = fc_TPP

# country summary data frame
c_sum_stats <- data.frame(GAUL = u_GAULS,
                          Country_burden_PY = NA,
                          Pop_density_release_area_median = NA,
                          Num_people_covered = NA,
                          Perc_national_pop_covered = NA,
                          Averted_med_costs_PY = NA,
                          Averted_VC_costs_PY = NA,
                          National_routine_VC_cost_PY = NA)

# decide on whether to return the minimum or median (cost to target all areas or median area cost)
sum_func_char <- "minimum" # choose "median" or "minimum"
if(sum_func_char == "median"){sum_func <- median}else{sum_func <- min}

# specific list of countries to extract detailed maps for:
f_countries <- data.frame(Name = c("Brazil", "Colombia", "Nigeria", "Bangladesh", "India", "Indonesia", "Sri Lanka"),
                          GAUL = c(37, 57, 182, 23, 115, 1013696, 231))
# focus country admin list
f_county_a2 <- list(ad2[ad2$COUNTRY_ID == "BRA", ],
                    ad2[ad2$COUNTRY_ID == "COL", ],
                    ad2[ad2$COUNTRY_ID == "NGA", ],
                    ad2[ad2$COUNTRY_ID == "BGD", ],
                    ad2[ad2$COUNTRY_ID == "IND", ],
                    ad2[ad2$COUNTRY_ID == "IDN", ],
                    ad2[ad2$COUNTRY_ID == "LKA", ])
names(f_county_a2) = f_countries$Name

# templates for country rasters
ID_vec = 1:length(pop_5km_v)
ID_vec[is.na(as.vector(baselayer))] = NA

# keep track of pixel IDs selected in each country to compile into a global raster
# also track their pricepoints
National_sel_IDs <- data.frame(ID = NA, GAUL = NA, Pop = NA, Averted_cases_mid = NA, Savings_5yr_mid = NA)
National_sel_IDs_50 <- data.frame(ID = NA, GAUL = NA, Pop = NA, Averted_cases_mid = NA, Savings_5yr_mid = NA)



# big extraction loop (loops through coutnry by country)

for(i in 1:length(u_GAULS)){
  # subset to specific country
  wp_cn <- world.pixels[world.pixels$GAUL == u_GAULS[i], ]
  # identify areas needed to meet national target
  cum_averted_cases_low = cumsum(wp_cn$Averted_cases_low)
  cum_averted_cases_mid = cumsum(wp_cn$Averted_cases_mid)
  cum_averted_cases_high = cumsum(wp_cn$Averted_cases_high)
  
  # find country burden match
  c_burd_match <- c_burd$Burd[match(u_GAULS[i], c_burd$GAUL)]
  
  pix_cut_low = which.max(cum_averted_cases_low >= (0.25 * c_burd_match))
  pix_cut_mid = which.max(cum_averted_cases_mid >= (0.25 * c_burd_match))
  pix_cut_high = which.max(cum_averted_cases_high >= (0.25 * c_burd_match))
  
  #if only half the reduction in global burden needs to be achieved using Wolbachia
  pix_cut50_low = which.max(cum_averted_cases_low >= (0.25 * 0.5 * c_burd_match))
  pix_cut50_mid = which.max(cum_averted_cases_mid >= (0.25 * 0.5 * c_burd_match))
  pix_cut50_high = which.max(cum_averted_cases_high >= (0.25 * 0.5 * c_burd_match))
  
  # save pixel IDs for global raster
  National_sel_IDs = rbind(National_sel_IDs, wp_cn[1:pix_cut_mid, c("ID", "GAUL", "Pop", "Averted_cases_mid", "Savings_5yr_mid")])
  National_sel_IDs_50 = rbind(National_sel_IDs_50, wp_cn[1:pix_cut50_mid, c("ID", "GAUL", "Pop", "Averted_cases_mid", "Savings_5yr_mid")])
  
  # calculate long term direct medical costs (i.e. programme break even cost) / people protected
  # to give target cost per person
  wp_cn$TPP_10yr_low = wp_cn$Savings_10yr_low / wp_cn$Pop
  wp_cn$TPP_10yr_mid = wp_cn$Savings_10yr_mid / wp_cn$Pop
  wp_cn$TPP_10yr_high = wp_cn$Savings_10yr_high / wp_cn$Pop
  
  wp_cn$TPP_5yr_low = wp_cn$Savings_5yr_low / wp_cn$Pop
  wp_cn$TPP_5yr_mid = wp_cn$Savings_5yr_mid / wp_cn$Pop
  wp_cn$TPP_5yr_high = wp_cn$Savings_5yr_high / wp_cn$Pop
  
  wp_cn$TPP_3yr_low = wp_cn$Savings_3yr_low / wp_cn$Pop
  wp_cn$TPP_3yr_mid = wp_cn$Savings_3yr_mid / wp_cn$Pop
  wp_cn$TPP_3yr_high = wp_cn$Savings_3yr_high / wp_cn$Pop
  
  # create country- specific TPP tables
  fc_TPP[i, 2:10] <-c(sum_func(wp_cn$TPP_10yr_low[1:pix_cut_low]),
                      sum_func(wp_cn$TPP_10yr_mid[1:pix_cut_mid]),
                      sum_func(wp_cn$TPP_10yr_high[1:pix_cut_high]),
                      sum_func(wp_cn$TPP_5yr_low[1:pix_cut_low]),
                      sum_func(wp_cn$TPP_5yr_mid[1:pix_cut_mid]),
                      sum_func(wp_cn$TPP_5yr_high[1:pix_cut_high]),
                      sum_func(wp_cn$TPP_3yr_low[1:pix_cut_low]),
                      sum_func(wp_cn$TPP_3yr_mid[1:pix_cut_mid]),
                      sum_func(wp_cn$TPP_3yr_high[1:pix_cut_high]))
  
  hc_TPP[i, 2:10] <-c(sum_func(wp_cn$TPP_10yr_low[1:pix_cut50_low]),
                      sum_func(wp_cn$TPP_10yr_mid[1:pix_cut50_mid]),
                      sum_func(wp_cn$TPP_10yr_high[1:pix_cut50_high]),
                      sum_func(wp_cn$TPP_5yr_low[1:pix_cut50_low]),
                      sum_func(wp_cn$TPP_5yr_mid[1:pix_cut50_mid]),
                      sum_func(wp_cn$TPP_5yr_high[1:pix_cut50_high]),
                      sum_func(wp_cn$TPP_3yr_low[1:pix_cut50_low]),
                      sum_func(wp_cn$TPP_3yr_mid[1:pix_cut50_mid]),
                      sum_func(wp_cn$TPP_3yr_high[1:pix_cut50_high]))
  
  # fill in country summary stats
  c_sum_stats$Country_burden_PY[i] = c_burd_match
  c_sum_stats$Num_people_covered[i] = sum(wp_cn$Pop[1:pix_cut_mid])
  c_sum_stats$Pop_density_release_area_median[i] = c_sum_stats$Num_people_covered[i] / (25 * pix_cut_mid)
  c_sum_stats$Perc_national_pop_covered[i] = 100 * sum(wp_cn$Pop[1:pix_cut_mid]) / sum(wp_cn$Pop)
  c_sum_stats$Averted_med_costs_PY[i] = sum(wp_cn$Savings_3yr_mid[1:pix_cut_mid] / 3)
  c_sum_stats$Averted_VC_costs_PY[i] = sum(wp_cn$Savings_Preferred_3yr_mid[1:pix_cut_mid] / 3)
  # need to deduct VC costs to calculate med costs only as "Savings_3yr_mid" includes both
  c_sum_stats$Averted_med_costs_PY[i] = c_sum_stats$Averted_med_costs_PY[i] - c_sum_stats$Averted_VC_costs_PY[i]
  c_sum_stats$National_routine_VC_cost_PY[i] = c_sum_stats$Averted_VC_costs_PY[i] / (0.35 * (3 / 12))
  
  # if in focus country list save a map
  if(u_GAULS[i] %in% f_countries$GAUL){
    # logical vector of whether pixels included in the release area
    ID_vec_S = ID_vec %in% wp_cn$ID[1:pix_cut_mid]
    ID_vec_S_50 = ID_vec %in% wp_cn$ID[1:pix_cut50_mid]
    
    # identifying the population, pre-intervention burden and costs supported in those release pixels
    ID_vec_S2 = match(ID_vec, wp_cn$ID[1:pix_cut_mid])
    ID_vec_S2A = wp_cn$TPP_5yr_mid[ID_vec_S2]
    ID_vec_S2B = wp_cn$Pop[ID_vec_S2]
    ID_vec_S2C = wp_cn$Averted_cases_mid[ID_vec_S2] / ((1 - Effectiveness))
    ID_vec_S2A[is.na(ID_vec_S2A)] = 0
    ID_vec_S2B[is.na(ID_vec_S2B)] = 0
    ID_vec_S2C[is.na(ID_vec_S2C)] = 0
    
    ID_vec_S2_50 = match(ID_vec, wp_cn$ID[1:pix_cut50_mid])
    ID_vec_S2_50A = wp_cn$TPP_5yr_mid[ID_vec_S2_50]
    ID_vec_S2_50B = wp_cn$Pop[ID_vec_S2_50]
    ID_vec_S2_50C = wp_cn$Averted_cases_mid[ID_vec_S2_50] / ((1 - Effectiveness))
    ID_vec_S2_50A[is.na(ID_vec_S2_50A)] = 0
    ID_vec_S2_50B[is.na(ID_vec_S2_50B)] = 0
    ID_vec_S2_50C[is.na(ID_vec_S2_50C)] = 0
    
    # identifying 
    
    examras1 = examras2 = examras3 = examras4 = examras5 = examras6 = examras7 = examras8 = baselayer
    # treatment areas
    values(examras1) = (values(examras1) + ID_vec_S - 1)
    values(examras2) = (values(examras2) + ID_vec_S_50 - 1)
    # supported costs in treatment areas
    values(examras3) = (values(examras3) + ID_vec_S2A - 1)
    values(examras4) = (values(examras4) + ID_vec_S2_50A - 1)
    # population
    values(examras5) = (values(examras5) + ID_vec_S2B - 1)
    values(examras6) = (values(examras6) + ID_vec_S2_50B - 1)
    # pre-intervention burden
    values(examras7) = (values(examras7) + ID_vec_S2C - 1)
    values(examras8) = (values(examras8) + ID_vec_S2_50C - 1)
    # crop to country extent
    examras1 = mask(crop(examras1, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras2 = mask(crop(examras2, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras3 = mask(crop(examras3, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras4 = mask(crop(examras4, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras5 = mask(crop(examras5, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras6 = mask(crop(examras6, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras7 = mask(crop(examras7, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    examras8 = mask(crop(examras8, ad0[ad0$GAUL_CODE == u_GAULS[i], ]), ad0[ad0$GAUL_CODE == u_GAULS[i], ])
    
    # save
    #writeRaster(examras1, file = paste0("/Users/eideobra/Dropbox/07_Indonesia/14_TPP/Figures/Target_",
    #                                    f_countries$Name[f_countries$GAUL == u_GAULS[i]],
    #                                    ".tif"),
    #            overwrite = T)
    #writeRaster(examras3, file = paste0("/Users/eideobra/Dropbox/07_Indonesia/14_TPP/Figures/Target_",
    #                                    f_countries$Name[f_countries$GAUL == u_GAULS[i]],
    #                                    "12_5perc.tif"),
    #            overwrite = T)
    
    # admin2 extract
    fc = (1:nrow(f_countries))[u_GAULS[i] == f_countries$GAUL]
    ad2_f <- exact_extract(examras1, f_county_a2[[fc]])
    ad2_f_50 <- exact_extract(examras2, f_county_a2[[fc]])
    ad2_f_cost <- exact_extract(examras3, f_county_a2[[fc]])
    ad2_f_cost_50 <- exact_extract(examras4, f_county_a2[[fc]])
    ad2_f_pop <- exact_extract(examras5, f_county_a2[[fc]])
    ad2_f_pop_50 <- exact_extract(examras6, f_county_a2[[fc]])
    ad2_f_burd <- exact_extract(examras7, f_county_a2[[fc]])
    ad2_f_burd_50 <- exact_extract(examras8, f_county_a2[[fc]])
    
    # aggregate extracted pixel results by municiplaity
    # more complex for cost because we need to choose the lowest supporting price pixel
    ad2_f <- lapply(ad2_f, function(x) c(sum(x[, 1] * x[, 2], na.rm = T), sum(x[, 2])))
    ad2_f_50 <- lapply(ad2_f_50, function(x) c(sum(x[, 1] * x[, 2], na.rm = T), sum(x[, 2])))
    ad2_f_cost_v <- rep(NA, length(ad2_f_cost))
    ad2_f_cost_v_50 <- rep(NA, length(ad2_f_cost_50))
    
    for(k in 1:length(ad2_f_cost)){
      temp = ad2_f_cost[[k]][, 1]
      temp[is.na(temp)] = 0
      if(any(temp > 0)){
        ad2_f_cost_v[k] = min(temp[temp > 0])
      }
      
      temp = ad2_f_cost_50[[k]][, 1]
      temp[is.na(temp)] = 0
      if(any(temp > 0)){
        ad2_f_cost_v_50[k] = min(temp[temp > 0])
      }
    }
    
    ad2_f_pop_v <- unlist(lapply(ad2_f_pop, function(x) sum(x[, 1] * x[, 2], na.rm = T)))
    ad2_f_pop_v_50 <- unlist(lapply(ad2_f_pop_50, function(x) sum(x[, 1] * x[, 2], na.rm = T)))
    ad2_f_burd_v <- unlist(lapply(ad2_f_burd, function(x) sum(x[, 1] * x[, 2], na.rm = T)))
    ad2_f_burd_v_50 <- unlist(lapply(ad2_f_burd_50, function(x) sum(x[, 1] * x[, 2], na.rm = T)))
    
    
    # reformat lists into matrices for area variables
    ad2_f = matrix(unlist(ad2_f), ncol = 2, byrow = T)
    ad2_f_50 = matrix(unlist(ad2_f_50), ncol = 2, byrow = T)
    
    ad2_tab <- data.frame(as.data.frame(f_county_a2[[fc]])[, 1:3],
                          Area_covered_KM2 = ad2_f[, 1] * 25,
                          Area_total_KM2 = ad2_f[, 2] * 25,
                          Cost_target = ad2_f_cost_v,
                          Population_covered = ad2_f_pop_v,
                          Burden_covered = ad2_f_burd_v)
    ad2_tab_50 <- data.frame(as.data.frame(f_county_a2[[fc]])[, 1:3],
                             Area_covered_KM2 = ad2_f_50[, 1] * 25,
                             Area_total_KM2 = ad2_f_50[, 2] * 25,
                             Cost_target = ad2_f_cost_v_50,
                             Population_covered = ad2_f_pop_v_50,
                             Burden_covered = ad2_f_burd_v_50)
    # trim to just municipalities selected
    ad2_tab = ad2_tab[ad2_tab$Area_covered_KM2 > 0, ]
    ad2_tab_50 = ad2_tab_50[ad2_tab_50$Area_covered_KM2 > 0, ]
    # save
    write.csv(ad2_tab, file = paste0("TPP_outputs/Munsum_",
                                     f_countries$Name[f_countries$GAUL == u_GAULS[i]],
                                     ".csv"))
    write.csv(ad2_tab_50, file = paste0("TPP_outputs/Munsum_",
                                        f_countries$Name[f_countries$GAUL == u_GAULS[i]],
                                        "12_5perc.csv"))
  }
}





### 10 summary and visualisaion of national level targets ###

# order countries by cost efficiency
fc_TPP = fc_TPP[order(fc_TPP$Y10_m, decreasing = T), ]
hc_TPP = hc_TPP[order(hc_TPP$Y10_m, decreasing = T), ]

# derive national TPPs
# defined here as price Wolbachia would need to be to reach national targets in 95% of dengue 
# endemic (> 10k symptomatic cases a year) countries

TPP_table <- rbind(fc_TPP[round(nrow(fc_TPP) *0.95, 0), ],
                   hc_TPP[round(nrow(hc_TPP) *0.95, 0), ])

row.names(TPP_table) <- c("Full contribution", "50% contribution")
TPP_table 
write.csv(TPP_table, file = "TPP_outputs/TPP_table_national.csv")

# trim global pixel lists to included countries, assign to rasters and save
National_sel_IDs = National_sel_IDs[National_sel_IDs$GAUL %in% fc_TPP[1:round(nrow(fc_TPP) *0.95, 0), "GAUL"], ]
National_sel_IDs_50 = National_sel_IDs_50[National_sel_IDs_50$GAUL %in% hc_TPP[1:round(nrow(hc_TPP) *0.95, 0), "GAUL"], ]

baselayer = ad0_ras >= 0
ID_vec = 1:length(pop_5km_v)
ID_vec[is.na(as.vector(baselayer))] = NA
ID_vec_100 = ID_vec %in% National_sel_IDs$ID
ID_vec_50 = ID_vec %in% National_sel_IDs_50$ID
cost_vec_100 = (National_sel_IDs$Savings_5yr_mid / National_sel_IDs$Pop)[match(ID_vec, National_sel_IDs$ID)]
cost_vec_100[is.na(cost_vec_100)] = 999999
cost_vec_50 = (National_sel_IDs_50$Savings_5yr_mid / National_sel_IDs_50$Pop)[match(ID_vec, National_sel_IDs_50$ID)]
cost_vec_50[is.na(cost_vec_50)] = 999999

examras1 = examras2 = examras3 = examras4 = baselayer
values(examras1) = (values(examras1) + ID_vec_100)
values(examras2) = (values(examras2) + ID_vec_50)
values(examras3) = cost_vec_100
values(examras4) = cost_vec_50

# save National targetting rasters
writeRaster(examras1, file = "TPP_outputs/Treat_areas_National_100perc_raster.tiff")
writeRaster(examras2, file = "TPP_outputs/Treat_areas_National_50perc_raster.tiff")

# processing municipality-level stats for country 12.5% and 25% goals
examras1 = examras1 == 2
examras2 = examras2 == 2
prop_covered_nat <- exact_extract(examras1, ad2, fun = "sum") / exact_extract(examras1^0, ad2, fun = "sum")
prop_covered_nat_12_5 <- exact_extract(examras2, ad2, fun = "sum") / exact_extract(examras2^0, ad2, fun = "sum")

# variable coding the same as other mun_sum files
mun_sum_national <- data.frame(A2_NAME = ad2$NAME,
                               A2_GAUL_CODE = ad2$GAUL_CODE,
                               A2COUNTRY = ad2$COUNTRY_ID,
                               Area_mun = mun_area,
                               Area_covered = mun_area * prop_covered_nat,
                               Pop_mun = exact_extract(pop_5km, ad2, fun = "sum"),
                               Pop_covered = exact_extract(pop_5km * examras1, ad2, fun = "sum"),
                               Burden_mun = exact_extract(denBurden, ad2, fun = "sum"),
                               Burden_covered = exact_extract(denBurden * examras1, ad2, fun = "sum"),
                               Cost_target = exact_extract(examras3, ad2, fun = "min"))
mun_sum_national = mun_sum_national[mun_sum_national$Area_covered > 0, ]

mun_sum_national_12_5 <- data.frame(A2_NAME = ad2$NAME,
                                    A2_GAUL_CODE = ad2$GAUL_CODE,
                                    A2COUNTRY = ad2$COUNTRY_ID,
                                    Area_mun = mun_area,
                                    Area_covered = mun_area * prop_covered_nat_12_5,
                                    Pop_mun = exact_extract(pop_5km, ad2, fun = "sum"),
                                    Pop_covered = exact_extract(pop_5km * examras2, ad2, fun = "sum"),
                                    Burden_mun = exact_extract(denBurden, ad2, fun = "sum"),
                                    Burden_covered = exact_extract(denBurden * examras2, ad2, fun = "sum"),
                                    Cost_target = exact_extract(examras4, ad2, fun = "min"))
mun_sum_national_12_5 = mun_sum_national_12_5[mun_sum_national_12_5$Area_covered > 0, ]
# and save
write.csv(mun_sum_national, file = "TPP_outputs/Munsum_national.csv")
write.csv(mun_sum_national_12_5, file = "TPP_outputs/Munsum_national_12_5perc.csv")


# processing country-level statistics
c_sum_stats = data.frame(Name = ad0$NAME[match(c_sum_stats$GAUL, ad0$GAUL_CODE)],
                         c_sum_stats)

# add in target implementation costs based on minimum and preferred global targets
# and national targets
c_sum_stats$Wol_Imp_cost_TPPmin = c_sum_stats$Num_people_covered * 2.33
c_sum_stats$Wol_Imp_cost_TPPpreferred = c_sum_stats$Num_people_covered * 0.24
c_sum_stats$National_Wol_target_min_PP = fc_TPP$Y5_m[match(c_sum_stats$GAUL, fc_TPP$GAUL)]
c_sum_stats$Wol_Imp_cost_TPP_National_min = c_sum_stats$National_Wol_target_min_PP * c_sum_stats$Num_people_covered

# order by cost efficiency
c_sum_stats = c_sum_stats[order(c_sum_stats$National_Wol_target_min_PP, 
                                decreasing = T), ]
# save
write.csv(c_sum_stats, file = "TPP_outputs/Country_characteristics.csv")






### 11 country cost figure ###

# distill country statistics
pf <- data.frame(Name= c_sum_stats$Name,
                 Price = c_sum_stats$National_Wol_target_min_PP,
                 Cum_burd = cumsum(c_sum_stats$Country_burden_PY),
                 Cum_averted_cases = cumsum(c_sum_stats$Country_burden_PY * 0.25))

# formatting
pf$Name = as.character(pf$Name)
pf$Name[pf$Name == "CÃ´te d'Ivoire"] = "Cote d'Ivoire"
c_r_tab <- read.csv("Country_region_table.csv")
pf$Region = c_r_tab$Region[match(pf$Name, c_r_tab$Name)]


# percentage increase in VC budget calculation
years_wol_cost_spread = 2
Bud_before = c_sum_stats$National_routine_VC_cost_PY * years_wol_cost_spread
# add emergency costs
Bud_before = Bud_before + Bud_before * 0.35 * (3 / 12)

Bud_after = Bud_before + 2.33  * c_sum_stats$Num_people_covered
Bud_after_5d = Bud_before + 5  * c_sum_stats$Num_people_covered
Bud_perc_change = 100* (Bud_after - Bud_before) / Bud_before
Bud_perc_change_5d = 100* (Bud_after_5d - Bud_before) / Bud_before
pf$Bud_perc_change = Bud_perc_change
pf$Bud_perc_change_5d = Bud_perc_change_5d

# names for just countries with major burdens or those that can support costs > 2.33
pf$Less_Name = as.character(pf$Name)
pf$Major_Name = as.character(pf$Name)
#pf$Less_Name[(c_sum_stats$National_Wol_target_min_PP < 2.33) & (c_sum_stats$Country_burden_PY < 500000)] = ""
pf$Less_Name[(c_sum_stats$National_Wol_target_min_PP < 2.33)] = ""
pf$Major_Name[(c_sum_stats$Country_burden_PY < 1000000)] = ""

# processing for linerange graphs
pf$CAC_L = c(0, pf$Cum_averted_cases[1:(nrow(pf) - 1)])
pf$CAC_M = pf$CAC_L + 0.5 * (pf$Cum_averted_cases - pf$CAC_L)
pf$Averted_cases = (pf$Cum_averted_cases - c(0, pf$Cum_averted_cases[1:(nrow(pf) - 1)])) / 1000000
# save data object to share
#write.csv(pf, file = "Figures/Country_prioritisation_data.csv")


p1A <- ggplot(pf, aes(x = Cum_averted_cases / 1000000, y = Price, label = Major_Name, group = Region)) +
  geom_linerange(aes(xmin = CAC_L / 1000000, xmax = Cum_averted_cases / 1000000, color = Region), size = 1) +
  geom_point(aes(color = Region)) +
  geom_text_repel(aes(color = Region), max.overlaps = 20) + 
  scale_y_log10(n.breaks = 10) +
  geom_hline(yintercept=2.33, lty = 2) +
  geom_hline(yintercept=0.44, lty = 2) +
  ylab("Averted medical and outbreak response costs per person (USD)") +
  xlab("Cumulative cases averted (millions)") + 
  theme_bw()

p1B <- ggplot(pf, aes(x = Cum_averted_cases / 1000000, y = Price, label = Less_Name, group = Region)) +
  geom_linerange(aes(xmin = CAC_L / 1000000, xmax = Cum_averted_cases / 1000000, color = Region), size = 1) +
  geom_point(aes(color = Region)) +
  geom_text_repel(aes(color = Region), max.overlaps = 20) + 
  scale_y_log10(n.breaks = 10) +
  geom_hline(yintercept=2.33, lty = 2) +
  geom_hline(yintercept=0.44, lty = 2) +
  ylab("Averted medical and outbreak response costs per person (USD)") +
  xlab("Cumulative cases averted (millions)") + 
  theme_bw()



# zoom in around the action zone 10-1$

p2 <- ggplot(pf, aes(x = Cum_averted_cases / 1000000, y = Price, label = Name, group = Region)) +
  geom_linerange(aes(xmin = CAC_L / 1000000, xmax = Cum_averted_cases / 1000000, color = Region), size = 1) +
  geom_point(aes(color = Region)) +
  geom_text_repel(aes(color = Region), max.overlaps = 20) + 
  geom_hline(yintercept=2.33, lty = 2) +
  geom_hline(yintercept=0.44, lty = 2) +
  scale_y_continuous(limits = c(0,3)) +
  scale_x_continuous(limits = c(14,25)) +
  ylab("Averted medical and outbreak response costs per person (USD)") +
  xlab("Cumulative cases averted (millions)") + 
  theme_bw()

p2 <- ggplot(pf, aes(x = Cum_averted_cases / 1000000, y = Price, label = Name, group = Region)) +
  geom_point(aes(color = Region)) +
  geom_text_repel(aes(color = Region)) + 
  geom_hline(yintercept=2.33, lty = 2) +
  geom_hline(yintercept=0.44, lty = 2) +
  scale_y_continuous(limits = c(0,3)) +
  scale_x_continuous(limits = c(14,25)) +
  ylab("Averted medical and outbreak response costs per person (USD)") +
  xlab("Cumulative cases averted (millions)") + 
  theme_bw()

# budget impact graphs
p3 <- ggplot(pf[pf$Price >= 2.33, ], aes(x = Bud_perc_change, y = Price, label = Major_Name, group = Region))+
  geom_point(aes(color = Region, size = Averted_cases)) +
  geom_text_repel(aes(color = Region)) + 
  scale_y_log10(n.breaks = 10) +
  ylab("Cost per person (USD) at which Wolbachia is cost neutral") +
  xlab("Percentage increase in VC budget to implement Wolbachia") + 
  theme_bw()

pf$Low_bud_Name = pf$Name
pf$Low_bud_Name[pf$Bud_perc_change > 100] = ""

p3B <- ggplot(pf[(pf$Price >= 2.33) & (pf$Bud_perc_change <= 100), ], 
              aes(x = Bud_perc_change, y = Price, label = Low_bud_Name, group = Region))+
  geom_point(aes(color = Region, size = Averted_cases)) +
  geom_text_repel(aes(color = Region), max.overlaps = 20) + 
  scale_y_log10(n.breaks = 10) +
  scale_x_continuous(limits = c(0,100)) +
  ylab("Cost per person (USD) at which Wolbachia is cost neutral") +
  xlab("Percentage increase in VC budget to implement Wolbachia") + 
  theme_bw()

p3C <- ggplot(pf[pf$Price >= 2.33, ], aes(x = Bud_perc_change_5d, y = Price, label = Major_Name, group = Region))+
  geom_point(aes(color = Region, size = Averted_cases)) +
  geom_text_repel(aes(color = Region)) + 
  scale_y_log10(n.breaks = 10) +
  ylab("Cost per person (USD) at which Wolbachia is cost neutral") +
  xlab("Percentage increase in VC budget to implement Wolbachia") + 
  theme_bw()

p3D <- ggplot(pf[(pf$Price >= 2.33) & (pf$Bud_perc_change <= 100), ], 
              aes(x = Bud_perc_change_5d, y = Price, label = Low_bud_Name, group = Region))+
  geom_point(aes(color = Region, size = Averted_cases)) +
  geom_text_repel(aes(color = Region), max.overlaps = 20) + 
  scale_y_log10(n.breaks = 10) +
  scale_x_continuous(limits = c(0,100)) +
  ylab("Cost per person (USD) at which Wolbachia is cost neutral") +
  xlab("Percentage increase in VC budget to implement Wolbachia") + 
  theme_bw()



# save maps
ggsave(p1A, filename = "TPP_outputs/Country_comparison_High_burden_countries.jpeg", width = 10, height = 10)
ggsave(p1B, filename = "TPP_outputs/Country_comparison_Higher_price_countries.jpeg", width = 10, height = 10)
ggsave(p2, filename = "TPP_outputs/Country_comparison_Lower_price_countries.jpeg", width = 10, height = 10)
ggsave(p3, filename = "TPP_outputs/Country_comparison_Budget_analysis.jpeg", width = 10, height = 5)
ggsave(p3B, filename = "TPP_outputs/Country_comparison_Budget_analysis_affordable.jpeg", width = 10, height = 5)
ggsave(p3C, filename = "TPP_outputs/Country_comparison_Budget_analysis_5dollars.jpeg", width = 10, height = 5)
ggsave(p3D, filename = "TPP_outputs/Country_comparison_Budget_analysis_affordable_5dollars.jpeg", width = 10, height = 5)


### END ###

