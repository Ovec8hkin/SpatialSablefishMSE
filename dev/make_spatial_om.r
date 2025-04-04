sable_om <- readRDS(file.path(here::here(), "data", "sablefish_om_big.RDS"))
spatial_dem_params <- readRDS(file.path(here::here(), "data", "spatial_dem_params_big.RDS"))
load(file.path(here::here(), "data", "spatial_sablefish_inputs.rda"))

trawl_fishery_selex <- sable_om$dem_params$sel[1,,,,2]
trawl_fishery_selex_mat <- afscOM::generate_param_matrix(trawl_fishery_selex, dimension_names = list("time"=1:200, "age"=2:31, "sex"=c("F", "M"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl")), by=c("age", "sex"), include_fleet_dim = TRUE)
spatial_dem_params$sel[,,,,2] <- trawl_fishery_selex_mat[,,,,2,drop=FALSE]


trawl_survey_selex <- sable_om$dem_params$surv_sel[1,,,,2]
trawl_survey_selex_mat <- afscOM::generate_param_matrix(trawl_survey_selex, dimension_names = list("time"=1:200, "age"=2:31, "sex"=c("F", "M"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl")), by=c("age", "sex"), include_fleet_dim = TRUE)
spatial_dem_params$surv_sel[,,,,2] <- trawl_survey_selex_mat[,,,,2,drop=FALSE]

dimension_names <- list("time"=1:200, "age"=2:31, "sex"=c("F", "M"), "region"=c("BS", "AI", "WGOA", "CGOA", "EGOA"), "fleet"=c("Fixed", "Trawl"))

names(sable_om)

init_naa <- aperm(spatial_sablefish_inputs$naa[,1,,,drop=FALSE], perm=c(2, 3, 4, 1))
dimnames(init_naa) = c(list("time"=1), dimension_names[2:4])

model_dimensions <- afscOM::get_model_dimensions(spatial_dem_params$sel)
model_options <- afscOM::setup_model_options(model_dimensions)
model_options$do_recruits_move <- FALSE
model_options$obs_pars <- list(
    # longline fishery, trawl fishery, longline survey, trawl survey
    is_survey   = c(0, 0, 1, 1),  # is this a survey (1) or fishery (0)
    qs          = c(1, 1, 6.113, 0.85), # catchability coefficient (q) for surveys
    catch_cv    = c(0.1, 0.1, 0, 0), # CV on catch for catch observations
    rpn         = c(0, 0, 1, 1), # should RPNs be computed (yes=1, no=0)
    rpn_cv      = c(0, 0, 0.1, 0.1), # RPN CV
    rpw         = c(0, 0, 1, 1), # should RPWs be computed (yes=1, no=0)
    rpw_cv      = c(0, 0, 0.1, 0.1), # RPW CV
    acs         = c(1, 1, 1, 1), # should age compositions be computed (yes=1, no=0)
    ac_samps    = c(200, 100, 200, 100), # total sample size for age composition observations
    ac_as_integers  = c(TRUE, TRUE, TRUE, TRUE), # return age comps as integers (TRUE) or proportions (FALSE)
    acs_agg_sex     = c(FALSE, FALSE, FALSE, FALSE) # should age comps be aggregated by sex
)

# Figure out average fleet/regional apportionment
catch <- spatial_sablefish_inputs$catch
catch <- aperm(catch, c(2,3,1))

regional_catch <- apply(catch, c(1, 3), sum)
regional_catch_props <- t(apply(regional_catch, 1, \(x) x/sum(x)))
fleet_catch_props <- array(NA, dim=dim(catch))
for(r in 1:5){
    for(y in 1:62){
        spat_catch <- regional_catch[y,r]
        flt_catch <- catch[y,,r]
        fleet_catch_props[y,,r] <- flt_catch/spat_catch
    }
}
recent_fleet_props <- afscOM::subset_matrix(fleet_catch_props, 40:62, d=1, drop=FALSE)
recent_regional_props <- afscOM::subset_matrix(regional_catch_props, 40:62, d=1, drop=FALSE)
region_fleet_catch_props <- apply(regional_catch_props, 2, mean)*t(apply(recent_fleet_props, c(2, 3), mean))

model_options$fleet_apportionment <- aperm(array(
    aperm(array(region_fleet_catch_props, dim=c(5, 2, 1)), c(3, 2, 1)),
    dim=c(2, 5, 200)
), c(3, 1 ,2))



spatial_sablefish_om <- list(
    dem_params = spatial_dem_params,
    model_options = model_options,
    init_naa = init_naa
)

saveRDS(spatial_sablefish_om, file.path(here::here(), "data", "spatial_sablefish_om.RDS"))
