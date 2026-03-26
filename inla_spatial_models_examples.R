library(sf)
library(spdep)
library(INLA)
library(ggplot2)

setwd("C:/Users/marika_dago/Desktop/2026_inla_spatial_models_oviedo")

# -------------------------------------------------------- -
# 1. Load North Carolina data ----
# -------------------------------------------------------- -
nc = st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

# Create area index
nc$ID.area = 1:nrow(nc)

# Example response/exposure:
# SID74 = SIDS (Sudden Infant Death Syndrome) cases
# BIR74 = births in 1974-78 --> expected cases (population proxy)
nc$E = nc$BIR74
nc$y = nc$SID74

# -------------------------------------------------------- -
# 2. Build adjacency graph for areal models ----
# -------------------------------------------------------- -
nb = poly2nb(nc)
nb2INLA("nc.adj", nb)


# -------------------------------------------------------- -
# 3. Create example covariates ----
# -------------------------------------------------------- -
# The NC dataset does not include enviromental variables,
# so here we create demo covariates for teaching purposes.

set.seed(123)

# "pm25": synthetic environmental exposure
nc$pm25 = abs(as.numeric(scale(nc$y)) + rnorm(nrow(nc), 0, 2))

# "temp": temperature
nc$temp = abs(as.numeric(scale(nc$y)) + rnorm(nrow(nc), 25, 4))


# -------------------------------------------------------- -
# 4. Fit spatial models ----
# -------------------------------------------------------- -

# 4a. Baseline model (only covariates) ----

m_base = inla(
  y ~ 1 + temp + pm25,
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_base)

# residual spatial pattern

# Pearson residuals
nc$fitted_base = m_base$summary.fitted.values$mean * nc$E
nc$pearson_base = (nc$y - nc$fitted_base) / sqrt(nc$fitted_base)

# Moran’s I test on residuals
lw = nb2listw(nb, style = "W")
moran.test(nc$pearson_base, lw)

# 4b. Random-effect model with covariates ----

# iid ----
m_iid = inla(
  y ~ 1 + temp + pm25 + f(ID.area, model = "iid"),
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_iid)

# 4c. Spatial models with covariates ----

# besag ----
m_besag = inla(
  y ~ 1 + temp + pm25 + f(ID.area, model = "besag", graph = nb),
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_besag)

nc$fitted_besag = m_besag$summary.fitted.values$mean * nc$E
nc$pearson_besag = (nc$y - nc$fitted_besag) / sqrt(nc$fitted_besag)
moran.test(nc$pearson_besag, lw)

# BYM ----
m_bym = inla(
  y ~ 1 + temp + pm25 + f(ID.area, model = "bym", graph = nb),
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_bym)

# BYM2 ----
m_bym2 = inla(
  y ~ 1 + temp + pm25 + f(ID.area, model = "bym2", 
                          graph = nb, scale.model = TRUE),
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_bym2)

# -------------------------------------------------------- -
# 5. Fit SVC models ----
# -------------------------------------------------------- -

# SVC besag ----
nc$ID.area2 = nc$ID.area

m_svc_besag = inla(
  y ~ 1 + temp + pm25 
    + f(ID.area, temp, model = "besag", graph = nb) +
    + f(ID.area2, pm25, model = "besag", graph = nb),
  family = "poisson",
  data = nc,
  E = E,
  control.compute = list(dic = TRUE, waic = TRUE),
  control.predictor = list(compute = TRUE)
)

summary(m_svc_besag)

nc$dev_temp = m_svc_besag$summary.random$ID.area$mean
nc$dev_pm25 = m_svc_besag$summary.random$ID.area2$mean


# -------------------------------------------------------- -
# 6. Models comparison ----
# -------------------------------------------------------- -

# Fit metrics ----
comp = data.frame(
  Model = c("base","iid", "Besag", "BYM", "BYM2", "SVC"),
  DIC   = c(m_base$dic$dic, m_iid$dic$dic, m_besag$dic$dic, 
            m_bym$dic$dic, m_bym2$dic$dic, m_svc_besag$dic$dic),
  WAIC  = c(m_base$waic$waic, m_iid$waic$waic, m_besag$waic$waic, 
            m_bym$waic$waic, m_bym2$waic$waicm, m_svc_besag$waic$waic)
)

comp

# Relative Risks ----
nc$rr_base   = exp(m_base$summary.linear.predictor$mean)
nc$rr_iid   = exp(m_iid$summary.linear.predictor$mean)
nc$rr_besag = exp(m_besag$summary.linear.predictor$mean)
nc$rr_bym   = exp(m_bym$summary.linear.predictor$mean)
nc$rr_bym2   = exp(m_bym2$summary.linear.predictor$mean)
nc$rr_svc   = exp(m_svc_besag$summary.linear.predictor$mean)


# --------------------- -
#### Maps ####
# --------------------- - 

# y
p = ggplot(nc) +
  geom_sf(aes(fill = y), color = "white", size = 0.2) +
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(fill = "SIDS")+
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(size = 12, face = "bold", color = "black"),
    strip.background = element_rect(fill = "grey80", color = "black", size = 0.5),
    legend.key.width = unit(0.5, "cm"),
    plot.title = element_text(size = 10)
  ) 

png(filename="y.png",
    width = 3000, height = 2600, units = "px", res=600)
p
dev.off()

# SVC spatial deviation temperature
p = ggplot(nc) +
  geom_sf(aes(fill = dev_temp), color = "grey80") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  theme_void() +
  labs(fill = expression(u[1]),
       title = "Temperature")
png(filename="dev_temp.png",
    width = 3000, height = 2600, units = "px", res=600)
p
dev.off()

# SVC spatial deviation PM2.5
p = ggplot(nc) +
  geom_sf(aes(fill = dev_pm25), color = "grey80") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  theme_void() +
  labs(fill = expression(u[2]),
       title = expression(PM[2.5]))
png(filename="dev_pm25.png",
    width = 3000, height = 2600, units = "px", res=600)
p
dev.off()

# MAPS RR

rr_vars = c("rr_base", "rr_iid", "rr_besag", "rr_bym", "rr_bym2", "rr_svc")
rr_labels = c("RR - Base", "RR - IID", "RR - Besag", "RR - BYM", "RR - BYM2", "RR - SVC")
rr_files = c("rr-base.png", "rr-iid.png", "rr-besag.png", "rr-bym.png", "rr-bym2.png", "rr-svc.png")

# common max for comparable scales across maps
rr_max = max(nc$rr_base, nc$rr_iid, nc$rr_besag, nc$rr_bym, nc$rr_bym2, nc$rr_svc, na.rm = TRUE)

# plotting function
make_rr_plot = function(data, fill_var, fill_label, rr_max) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_var]]), 
            color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "plasma", limits = c(0, rr_max)) +
    labs(fill = fill_label) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", 
                                          color = NA),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "white", 
                                         color = NA),
          strip.text = element_text(size = 12, 
                                    face = "bold", 
                                    color = "black"),
          strip.background = element_rect(fill = "grey80", 
                                          color = "black", 
                                          linewidth = 0.5),
          legend.key.width = unit(0.5, "cm"),
          plot.title = element_text(size = 10))
}

# loop 
for (i in seq_along(rr_vars)) {
  p = make_rr_plot(data = nc,
                   fill_var = rr_vars[i],
                   fill_label = rr_labels[i],
                   rr_max = rr_max)
  
  png(filename = rr_files[i], 
      width = 3000, height = 2600, units = "px", res = 600)
  print(p)
  dev.off()
}
