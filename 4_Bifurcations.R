source("0_Functions.R")

# limit cycle

ch <- 0.5#3
chm <- 10#3 # needs to be bigger than or equal to cd
dh <- 1
dhp <- 0#1.25

cp <- 5
dp <- 1

cm <- 10#4
dm <- 1

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

ini_state <- c(1, 0.5, 0.5) * 0.125
names(ini_state) <- c("Host", "Pathogen", "Mutualist")

end_time <- 50
time_step <- 0.001

out_dyn <- IntegrateDynamics(ini_state, pars,
                             end_time, time_step,
                             fn = CoralMutPathDynamics)

melt_dyn <- melt(out_dyn, id.vars = c("time"))

plDyn <- ggplot(melt_dyn, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_classic() + scale_y_log10() +
  labs(x = "Time", y = "Frequency", color = "", linetype = "") +
  ggtitle("A") +
  scale_color_manual(values = c("darkblue", "darkred", "darkgreen")) +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 10),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plDyn

### bifurcation

ch <- 0.5
chm <- c(3.56, 3.57) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

end_time <- 2e3
time_step <- 0.1

out_data <- data.frame()

for(i in 1:nrow(iterated_params)) {
  
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)
  
  root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
  h_soln <- root_soln$root
  mp_soln <- GetPM(h_soln, cur_pars)
  p_soln <- mp_soln[1]
  m_soln <- mp_soln[2]
  
  ini_state <- c(h_soln, p_soln, m_soln) + rnorm(n = 3, mean = 0, sd = 0.001)
  ini_state[ini_state <= 0] <- 0.001
  names(ini_state) <- c("host", "pathogen", "mutualist")
  
  out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                               end_time, time_step,
                               fn = CoralMutPathDynamics)
  
  cur_dyn <- cbind(cur_params, out_dyn)
  
  out_data <- rbind(out_data, cur_dyn)
  
}

end_points <- out_data %>%
  filter(time == max(time))

plot_data <- out_data %>%
  mutate(PatEnd = case_when(chm == min(chm) ~ end_points[end_points$chm == min(chm),]$pathogen,
                            chm == max(chm) ~ end_points[end_points$chm == max(chm),]$pathogen)) %>%
  mutate(MutEnd = case_when(chm == min(chm) ~ end_points[end_points$chm == min(chm),]$mutualist,
                            chm == max(chm) ~ end_points[end_points$chm == max(chm),]$mutualist))

plBifurcation <- ggplot(plot_data, aes(x = pathogen, y = mutualist, color = time)) +
  geom_path(linewidth = 0.5) +
  geom_point(aes(x = PatEnd, y = MutEnd), color = "red", size = 2) +
  facet_wrap(~chm, scales = "free", labeller = label_bquote(cols = c[hm] == .(chm)), nrow = 1) +
  theme_classic() +
  scale_color_gradient(low = "darkblue", high = "green") +
  labs(x = "Pathogen",
       y = "Mutualist",
       color = "Time") +
  ggtitle("B") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plBifurcation

### allee effect

ch <- 0.5
chm <- 5 # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- 10
dm <- 1

in_pars <- data.frame(ch = ch, chm = chm, dh = dh, dhp = dhp,
                      cp = cp, dp = dp,
                      cm = cm, dm = dm)

end_time <- 5e2
time_step <- 1
freq_cutoff <- 1e-7

ini_hosts <- seq(0.001, 0.9, length.out = 100)
ini_muts <- seq(0.0001, 0.15, length.out = 100)

ini_path <- 0.0001 # should be less than the smaller value above

out_data <- data.frame()

for(i in 1:nrow(in_pars)) {
  
  for(ini_host in ini_hosts) {
    for(ini_mut in ini_muts) {
      
      cur_params <- in_pars[i,]
      cur_pars <- as.list(cur_params)
      
      ini_state <- c(ini_host, ini_path, ini_mut)
      names(ini_state) <- c("Host", "Pathogen", "Mutualist")
      
      out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                                   end_time, time_step,
                                   fn = CoralMutPathDynamics)
      
      end_state <- as.numeric(out_dyn[nrow(out_dyn),2:4])
      
      cur_params$IniMut <- ini_mut
      cur_params$IniHost <- ini_host
      cur_params$IniPath <- ini_path
      
      out_coexist <- data.frame(Coexist = ifelse(prod(end_state > freq_cutoff),
                                                 "Coexistence",
                                                 "No coexistence"))
      
      cur_dyn <- cbind(cur_params, out_coexist)
      out_data <- rbind(out_data, cur_dyn)
    }
  }
  
}

out_data <- out_data %>%
  mutate(Coexist = ifelse(IniMut + IniPath > IniHost, NA, Coexist))

plIniCond <- ggplot(out_data,
                    aes(x = IniHost, y = IniMut, fill = Coexist)) +
  geom_tile() + theme_classic() +
  scale_fill_manual(values = c("blue", "gray"), na.value="white", na.translate = F) +
  labs(x = "Initial Frequency of the Host",
       y = "Initial Frequency of the Mutualist",
       fill = "") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("C") +
  theme(axis.text = element_text( size = 10 ),
        legend.position = "top",
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plIniCond

layout_vec <- c(1, 1, 2, 2, 2, 3, 3)
jpeg("./figs/Fig4Bifurcations.jpeg",
     width = 5000, height = 1250, res = 300)
grid.arrange(plDyn, plBifurcation, plIniCond, nrow = 1, layout_matrix = matrix(layout_vec, nrow = 1, ncol = length(layout_vec)))
dev.off()




