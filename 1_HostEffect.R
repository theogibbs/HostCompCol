source("0_Functions.R")

# Feasibility and stability of coexistence state as host mortality and colonization rates are varied ----
ch <- seq(1.5, 10, length.out = 200)
dh <- seq(0.0001, 1.5, length.out = 200)
dhp <- 0

cp <- 5
dp <- 1

cm <- c(6, 10, 25)
dm <- 1

in_pars <- crossing(ch = ch, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)
in_pars$chm <- in_pars$ch

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()
plan(multisession)
out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind',
  .inorder = FALSE) %dofuture% {
  
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)
  
  if (cur_pars$ch > cur_pars$dh) {
  root_soln <- uniroot(f = PredEq, interval = c(1e-10, 1), cur_pars)
  h_soln <- root_soln$root
  mp_soln <- GetPM(h_soln, cur_pars)
  p_soln <- mp_soln[1]
  m_soln <- mp_soln[2]
  } else {
    h_soln = 0
    p_soln = 0
    m_soln = 0
  }

  J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
  cur_eig <- GetEig(J)
  
  cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
  cur_stable <- cur_eig < 0
  
  cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                        ifelse(cur_stable, "\nFeasible\nand stable\n",
                               "\nFeasible\nbut unstable\n"))
  
  cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
  cur_dyn
}

plHeatMapHost <- ggplot(out_data,
                         aes(x = dh, y = ch, fill = Outcome)) +
  facet_wrap(~cm, labeller = label_bquote(c[m] == .(cm))) +
  geom_tile() + theme_classic() +
  labs(x = expression("Host Mortality" ~ (d[h])),
       y = expression("Host Colonization" ~ (c[h])),
       fill = "") +
  scale_fill_manual("Coexistence\nstatus:", values = c("#0072B2", "#D55E00")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  ggtitle("C")
plHeatMapHost # Figure 1 C in the main text

# Example trajectories as host mortality is varied ----
ch <- 2
chm <- ch # needs to be bigger than or equal to ch
dh <- c(1, 1.1, 1.2)
dhp <- 0

cp <- 5
dp <- 1

cm <- 6
dm <- 1

in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

end_time <- 50
time_step <- 0.01

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
  select(c("dh", "time", "Host", "Pathogen", "Mutualist")) %>%
  melt(id.vars = c("time", "dh"))

plDyn <- ggplot(melt_dyn, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_classic() +
  facet_wrap(~dh, scales = "free", labeller = label_bquote(cols = d[h] == .(dh))) +
  labs(x = "Time", y = "Frequency", color = "") +
  ggtitle("A") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.05))+
  scale_color_viridis_d() +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plDyn # Figure 1 A in the main text

# Feasibility and stability of coexistence state as host mortality and colonization rates are varied ----
ch <- 2
dh <- c(0, 0.5, 1)
dhp <- 0

cp <- seq(5, 10, length.out = 200)
dp <- 1

cm <- seq(5, 35, length.out = 200)
dm <- 1

in_pars <- crossing(ch = ch, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)
in_pars$chm <- in_pars$ch

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

out_data <- data.frame()
plan(multisession)
out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind',
  .inorder = FALSE) %dofuture% {
    
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)
    
    if (cur_pars$ch > cur_pars$dh) {
      root_soln <- uniroot(f = PredEq, interval = c(1e-10, 1), cur_pars)
      h_soln <- root_soln$root
      mp_soln <- GetPM(h_soln, cur_pars)
      p_soln <- mp_soln[1]
      m_soln <- mp_soln[2]
    } else {
      h_soln = 0
      p_soln = 0
      m_soln = 0
    }
    
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    cur_eig <- GetEig(J)
    
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    cur_stable <- cur_eig < 0
    
    cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                          ifelse(cur_stable, "\nFeasible\nand stable\n",
                                 "\nFeasible\nbut unstable\n"))
    
    cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
    cur_dyn
  }

plHeatMapClassicComp <- ggplot(out_data,
                        aes(x = cp, y = cm, fill = Outcome)) +
  facet_wrap(~dh, labeller = label_bquote(d[h] == .(dh))) +
  geom_tile() + theme_classic() +
  labs(x = expression("Pathogen Colonization" ~ (c[p])),
       y = expression("Mutualist Colonization" ~ (c[m])),
       fill = "") +
  scale_fill_manual("Coexistence\nstatus:", values = c("#0072B2", "#D55E00")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot") +
  ggtitle("B")
plHeatMapClassicComp # Figure 1 B in the main text


# Save figure ----
jpeg("./figs/Fig1HostEffect.jpeg",
     width = 3000, height = 3000, res = 300)
grid.arrange(plDyn, plHeatMapClassicComp, plHeatMapHost, nrow = 3)
dev.off()



