source("0_Functions.R")

ch <- c(2, 3, 4)
chm <-  -1 # needs to be bigger than or equal to ch
dh <- 1
dhp <- seq(0, 2, length.out = 100)

cp <- 5
dp <- 1

cm <- seq(2.5, 8, length.out = 100)
dm <- 1

in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)
in_pars$chm <- in_pars$ch

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()

for(i in 1:nrow(iterated_params)) {
  
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)
  
  root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
  h_soln <- root_soln$root
  mp_soln <- GetPM(h_soln, cur_pars)
  p_soln <- mp_soln[1]
  m_soln <- mp_soln[2]
  
  J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
  cur_eig <- GetEig(J)
  
  cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
  cur_stable <- cur_eig < 0
  
  cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                        ifelse(cur_stable, "\nFeasible\nand stable\n",
                               "\nFeasible\nbut unstable\n"))
  
  cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
  out_data <- rbind(out_data, cur_dyn)
  
}

plHeatMapDeath <- ggplot(out_data,
       aes(x = dhp, y = cm, fill = Outcome)) +
  facet_wrap(~ch, labeller = label_bquote(c[h] == .(ch))) +
  geom_tile() + theme_classic() +
  labs(x = expression("Added Host Mortality from the Pathogen" ~ (d[hp])),
       y = expression("Mutualist Colonization"~(c[m])),
       fill = "") +
  scale_fill_manual(values = c("green", "red")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  ggtitle("B")
plHeatMapDeath


ch <- 4
chm <- 4 # needs to be bigger than or equal to ch
dh <- 1
dhp <- c(0, 1.5, 2)

cp <- 5
dp <- 1

cm <- 4
dm <- 1

in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

end_time <- 50
time_step <- 0.001

dyn_data <- data.frame()

for(i in 1:nrow(iterated_params)) {
  
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)
  
  ini_state <- c(0.2, 0.05, 0.05)
  names(ini_state) <- c("Host", "Pathogen", "Mutualist")
  
  out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                               end_time, time_step,
                               fn = CoralMutPathDynamics)
  
  cur_dyn <- cbind(cur_params, out_dyn)
  
  dyn_data <- rbind(dyn_data, cur_dyn)
  
}

melt_dyn <- dyn_data %>%
  select(c("dhp", "time", "Host", "Pathogen", "Mutualist")) %>%
  melt(id.vars = c("time", "dhp"))

plDyn <- ggplot(melt_dyn, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_classic() + scale_y_log10() +
  facet_wrap(~dhp, scales = "free", labeller = label_bquote(cols = d[hp] == .(dhp))) +
  labs(x = "Time", y = "Frequency", color = "") +
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

jpeg("./figs/Fig2PathDeath.jpeg",
     width = 3000, height = 2000, res = 300)
grid.arrange(plDyn, plHeatMapDeath, nrow = 2)
dev.off()
