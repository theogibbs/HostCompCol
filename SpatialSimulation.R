source("0_Functions.R")

set.seed(1)

ch <- 1.5
chm <- 10  # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

N <- 5e2
Rh <- 10
Rb <- 1

host_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rh)
bact_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rb)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm,
             N = N, Rh = Rh, Rb = Rb,
             host_kern = host_kern, bact_kern = bact_kern)

ini_state <- sample(x = c("Empty",
                          "Host",
                          "Pathogen",
                          "Mutualist"
                          ),
                    size = N, replace = TRUE)
ini_state <- data.frame(Site = 1:pars$N, Occ = ini_state)

end_step <- 1e5

collect_length <- 1000
collect_interval <- end_step / collect_length
collect_spatial_steps <- seq(0, end_step, by = collect_interval)

collect_freq_steps <- 1:end_step

out_list <- RunSpatialSim(pars = pars,
                          ini_state = ini_state,
                          end_step = end_step,
                          collect_spatial_steps = collect_spatial_steps,
                          collect_freq_steps = collect_freq_steps,
                          print_time = TRUE)

out_freq <- out_list$freq
out_freq <- as.data.frame(out_freq)
colnames(out_freq) <- c("Step", "Time", "Host", "Pathogen", "Mutualist")

melt_freq <- out_freq %>%
  melt(id.vars = c("Step", "Time"))
colnames(melt_freq) <- c("Step", "Time", "Occ", "Frequency")

melt_freq$Occ <- factor(melt_freq$Occ, levels = c("Host", "Pathogen", "Mutualist"))

plFreqDyn <- ggplot(melt_freq, aes(x = Time, y = Frequency, color = Occ)) +
  geom_line() + theme_classic() + scale_y_log10() +
  labs(x = "Time", y = "Frequency", color = "") +
  scale_color_manual(values = c("darkblue", "darkred", "darkgreen"))
plFreqDyn

out_data <- out_list$spatial
out_data$Occ <- factor(out_data$Occ, levels = c("Empty", "Host",
                                                "Pathogen", "Mutualist"))

out_data$Occ[out_data$Occ == "Empty"] <- NA

ini_host <- sum(ini_state$Occ == "Host")
ini_path <- sum(ini_state$Occ == "Pathogen")
ini_mut <- sum(ini_state$Occ == "Mutualist")

ini_state <- c(ini_host + ini_path + ini_mut, ini_path, ini_mut) / pars$N
names(ini_state) <- c("Host", "Pathogen", "Mutualist")

end_time <- max(out_freq$Time)
time_step <- end_time / 1000
out_dyn <- IntegrateDynamics(ini_state, pars,
                             end_time, time_step,
                             fn = CoralMutPathDynamics)

melt_dyn <- melt(out_dyn, id.vars = c("time"))

plot_freq <- melt_freq %>%
  select(!c("Step"))
plot_freq$Source <- "Spatial"

colnames(melt_dyn) <- c("Time", "Occ", "Frequency")
melt_dyn$Source <- "Deterministic"

melt_dyn <- rbind(melt_dyn, plot_freq) %>%
  mutate(Source = factor(Source, levels = c("Spatial", "Deterministic")))

plComp <- ggplot(melt_dyn, aes(x = Time, y = Frequency, color = Occ, linetype = Source)) +
  geom_line(linewidth = 1) + 
  theme_classic() + 
  labs(x = "Time", y = "Frequency", color = "", linetype = "") +
  ggtitle("B") +
  scale_color_viridis_d(guide = "none") +
  theme(text = element_text(size=20),
        legend.text=element_text(size = 15),
        legend.position = "top",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  ylim(c(0, 1))
plComp

diff_time <- c(0, unique(out_data$Time))
diff_time <- diff(diff_time)
diff_time <- rep(diff_time, each = pars$N)
out_data$DiffTime <- diff_time
out_data$DiffTime <- out_data$DiffTime * 2

plSpatialDynamics <- ggplot(out_data, aes(x = Time, y = Site,
                                          width = DiffTime, fill = Occ)) +
  geom_tile() + theme_classic() +
  scale_y_continuous(limits = c(0, N+1 + 5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.25, max(out_data$Time) + 0.25), expand = c(0, 0)) +
  scale_fill_viridis_d(na.translate = FALSE) +
  theme(#axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.position = "top",
    #axis.ticks.y = element_blank(),
    text = element_text(size=20),
    #axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.text=element_text(size = 15),
    plot.title.position = "plot",
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
  ggtitle("A") + labs(fill = "")
plSpatialDynamics

sim_freqs <- read_csv2(file = "../ManuscriptHostCompCol/simdata/sim0213_spatialcomp.csv")

melt_freqs <- sim_freqs %>%
  select(c("Step", "Time", "Host", "Pathogen", "Mutualist", "Spatial", "cm", "ch", "chm")) %>%
  melt(id.vars = c("Step", "Time", "Spatial", "cm", "ch", "chm")) %>%
  group_by(Spatial, cm, ch, chm, variable) %>%
  summarise(MeanVal = median(value),
            LowerVal = quantile(value, probs = 0.25),
            UpperVal = quantile(value, probs = 0.75)) %>%
  mutate(MutBenefit = ifelse(chm > ch, "c[hm] == 10", "c[hm] == c[h]")) %>%
  mutate(Spatial = ifelse(Spatial == "NON-SPATIAL", "Stochastic", "Spatial")) %>%
  mutate(Spatial = fct_rev(Spatial))


pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

ode_chs <- seq(min(sim_freqs$ch), max(sim_freqs$ch), length.out = 100)
cms <- unique(sim_freqs$cm)
in_chms <- 10
ode_preds <- data.frame()

for(cur_cm in cms) {
  for(cur_ch in ode_chs) {
    chms <- c(cur_ch, in_chms[in_chms >= cur_ch])
    for(cur_chm in chms) {
      
      pars$cm <- cur_cm
      pars$ch <- cur_ch
      pars$chm <- cur_chm
      
      # solve for equilibrium for each population
      
      if(cur_ch < cur_chm) {
        root_soln <- uniroot(f = PredEq, interval = c(0, 1), pars)
        h_soln <- root_soln$root
      } else {
        h_soln <- 1 - pars$dh / pars$ch
      }
      mp_soln <- GetPM(h_soln, pars)
      p_soln <- mp_soln[1]
      m_soln <- mp_soln[2]
      
      if(p_soln < 0) {
        m_soln <- h_soln - (pars$dh + pars$dm) / pars$cm
      }
      
      if(m_soln < 0) {
        p_soln <- h_soln - (pars$dh + pars$dp) / pars$cp
      }
      
      # determine feasibility
      cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
      
      J <- BuildJacobian(h_soln, p_soln, m_soln, pars)
      cur_eig <- GetEig(J)
      
      # label outcomes accordingly with respective population, feasibility, and stability.
      cur_preds <- data.frame(Host = h_soln, Pathogen = p_soln, Mutualist = m_soln,
                              cm = cur_cm, ch = cur_ch, chm = cur_chm,
                              Feasibility = cur_feas, Eigenvalue = cur_eig)
      
      ode_preds <- rbind(ode_preds, cur_preds)
    }
  }
}

melt_ode_preds <- ode_preds %>%
  select(!c("Feasibility", "Eigenvalue")) %>%
  melt(id.vars = c("cm", "ch", "chm")) %>%
  mutate(value = pmax(0, value)) %>%
  mutate(MutBenefit = ifelse(chm > ch, "c[hm] == 10", "c[hm] == c[h]")) %>%
  mutate(Spatial = "Deterministic")


melt_freqs$MutBenefit <- factor(melt_freqs$MutBenefit, levels = c("c[hm] == c[h]", "c[hm] == 10"))

plPercolation <- ggplot(melt_freqs, aes(x = ch, y = MeanVal, shape = Spatial, color = Spatial)) +
  geom_line(data = melt_ode_preds, aes(x = ch, y = value, color = Spatial)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerVal, ymax = UpperVal, width = 0.15)) +
  theme_classic() +
  scale_color_viridis_d(option = "plasma", end = 0.75) +
  facet_grid(MutBenefit~variable,
             labeller = label_parsed, 
             scales = "free") +
  theme(text = element_text(size=20),
        panel.spacing = unit(2.5, "lines"),
        legend.text=element_text(size = 15),
        legend.position = "top",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = expression("Host Colonization" ~ (c[h])),
       y = "Frequency",
       color = " ",
       shape = " ") +
  ggtitle("C")
plPercolation


layout_mat <- rbind(c(1, 3, 3), c(2, 3, 3))
# writing out graph
jpeg("./figs/FigSpace.jpeg",
     width = 4000, height = 2250, res = 300)
grid.arrange(plSpatialDynamics, plComp, plPercolation, layout_matrix = layout_mat)
dev.off()

