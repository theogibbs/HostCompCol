source("0_Functions.R")

ch <- 0
chm <- 0  # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

N <- 5e2
Rh <- 1
Rb <- 1

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm,
             N = N, Rh = Rh, Rb = Rb)

ini_state <- sample(x = c("Empty",
                          "Pathogen",
                          "Host",
                          "Mutualist"
),
size = pars$N, replace = TRUE)
ini_state <- data.frame(Site = 1:pars$N, Occ = ini_state)

end_step <- 1e5

collect_spatial_steps <- NA

collect_freq_steps <- (end_step - 100):end_step

#out_list <- RunSpatialSim(pars = pars,
#                          ini_state = ini_state,
#                          end_step = end_step,
#                          collect_spatial_steps = collect_spatial_steps,
#                          collect_freq_steps = collect_freq_steps)
#out_list$freq
spatial_settings <- c("SPATIAL", "NON-SPATIAL")
cms <- 10
chs <- seq(0.05, 2.25, length.out = 7)
in_chms <- 10

start_time <- Sys.time()

out_freqs <- data.frame()
num_repl <- 10

for(i in 1:num_repl) {
  print(paste("REPLICATE:", i, "OUT OF", num_repl, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"))
  for(cur_space in spatial_settings) {
    print(cur_space)
    if(cur_space == "NON-SPATIAL") {
      Rh <- pars$N
      Rb <- pars$N
      
      host_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rh)
      bact_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rb)
      
      pars$host_kern <- host_kern
      pars$bact_kern <- bact_kern
      
    } else {
      Rh <- 10
      Rb <- 1
      
      host_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rh)
      bact_kern <- BuildKernel("1D Lattice", N = N, kern_length = Rb)
      
      pars$host_kern <- host_kern
      pars$bact_kern <- bact_kern
      
    }
    for(cur_cm in cms) {
      print(paste("Current mutualist colonization:", round(cur_cm, 3)))
      for(cur_ch in chs) {
        print(paste("   Current host colonization:", round(cur_ch, 3)))
        chms <- c(cur_ch, in_chms[in_chms >= cur_ch])
        for(cur_chm in chms) {
          print(paste("       Current host colonization from the mutualist:", round(cur_chm, 3)))
          
          pars$cm <- cur_cm
          pars$ch <- cur_ch
          pars$chm <- cur_chm
          
          out_list <- RunSpatialSim(pars = pars,
                                    ini_state = ini_state,
                                    end_step = end_step,
                                    collect_spatial_steps = collect_spatial_steps,
                                    collect_freq_steps = collect_freq_steps,
                                    print_time = FALSE)
          
          cur_freq <- out_list$freq
          cur_freq$Spatial <- cur_space
          cur_freq$cm <- cur_cm
          cur_freq$ch <- cur_ch
          cur_freq$chm <- cur_chm
          cur_freq$Replicate <- i
          out_freqs <- rbind(out_freqs, cur_freq)
          
        }
      }
    }
  }
}


Sys.time() - start_time

melt_freqs <- out_freqs %>%
  select(!c("Replicate")) %>%
  melt(id.vars = c("Step", "Time", "Spatial", "cm", "ch", "chm")) %>%
  group_by(Spatial, cm, ch, chm, variable) %>%
  summarise(MeanVal = mean(value),
            LowerVal = quantile(value, probs = 0.1),
            UpperVal = quantile(value, probs = 0.9)) %>%
  mutate(MutBenefit = ifelse(chm > ch, "c[hm] == 10", "c[hm] == c[h]")) %>%
  mutate(Spatial = ifelse(Spatial == "NON-SPATIAL", "Stochastic", "Spatial")) %>%
  mutate(Spatial = fct_rev(Spatial))


pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

ode_chs <- seq(min(chs), max(chs), length.out = 100)
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


plPercolation <- ggplot(melt_freqs, aes(x = ch, y = MeanVal, color = Spatial)) +
  geom_line(data = melt_ode_preds, aes(x = ch, y = value, color = Spatial)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LowerVal, ymax = UpperVal, width = 0.2)) +
  theme_classic() +
  scale_color_viridis_d(option = "turbo") +
  facet_wrap(MutBenefit~variable,
             labeller = label_parsed, 
             scales = "free") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  labs(x = expression("Host Colonization" ~ (c[h])),
       y = "Frequency",
       color = " ")
plPercolation

# writing out graph
jpeg("./figs/FigSpatialComparison.jpeg",
     width = 2500, height = 1500, res = 300)
plPercolation
dev.off()

in_pars <- data.frame(dh = dh, dhp = dhp,
                      cp = dp, dp = dp,
                      dm = dm,
                      N = N, Rh = Rh, Rb = Rb)

write_data <- cbind(out_freqs, in_pars)

write_csv2(write_data, file = "../ManuscriptHostCompCol/simdata/sim0213_spatialcomp.csv")


