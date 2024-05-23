source("0_Functions.R")

ch <- 0.5
chm <- seq(0.5, 10, length.out = 200) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- c(10, 12)
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

ini_cond_sds <- c(0.001, 0.5)

end_time <- 500
time_step <- 1

out_data <- data.frame()

for(i in 1:nrow(iterated_params)) {
  for(ini_cond_sd in ini_cond_sds) {
    
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)
    
    root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
    h_soln <- root_soln$root
    mp_soln <- GetPM(h_soln, cur_pars)
    p_soln <- mp_soln[1]
    m_soln <- mp_soln[2]
    
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    
    ini_state <- c(h_soln, p_soln, m_soln) + rnorm(n = 3, mean = 0, sd = ini_cond_sd)
    ini_state[ini_state <= 0] <- 0.0001
    names(ini_state) <- c("Host", "Pathogen", "Mutualist")
    
    out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                                 end_time, time_step,
                                 fn = CoralMutPathDynamics)
    
    out_dyn <- out_dyn %>%
      filter(time > end_time / 2) %>%
      melt(id.vars = c("time"))
    
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    cur_eig <- GetEig(J)
    
    cur_params$IniCondSd <- ini_cond_sd
    
    cur_dyn <- cbind(cur_params, out_dyn) %>%
      mutate(pred_soln = case_when(variable == "Host" ~ h_soln,
                                   variable == "Mutualist" ~ m_soln,
                                   variable == "Pathogen" ~ p_soln)) %>%
      mutate(Eigenvalue = cur_eig,
             Stable = ifelse(cur_feas, ifelse(cur_eig < 0,
                                              "Feasible and stable",
                                              "Feasible but unstable"),
                             "Not feasible"))
    
    out_data <- rbind(out_data, cur_dyn)
  }
  
}

proc_data <- out_data %>%
  mutate(pred_soln = ifelse(pred_soln > 0, pred_soln, NA))

#ggplot(out_data,
#       aes(x = chm, y = Eigenvalue, color = as.factor(cm))) +
#  geom_point() + theme_classic() +
#  geom_hline(yintercept = 0) +
#  labs(x = "Added Host Colonization from the Mutualist",
#       y = "Eigenvalue",
#       color = "Mutualist\nColonization")

plMutBenefit <- ggplot(proc_data,
       aes(x = chm, y = value)) +
  facet_grid(cm~variable, labeller = label_bquote(rows = c[m] == .(cm), cols = .(variable))) +
  geom_point(size = 2, alpha = 0.005) + theme_classic() +
  geom_line(linewidth = 1, aes(x = chm, y = pred_soln, color = Stable)) +
  scale_color_manual(values = c("green", "orange", "red")) +
  labs(x = expression("Added Host Colonization from the Mutualist" ~ (c[hm])),
       y = "Frequency",
       color = "",
       linetype = "") +
  ggtitle("A") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        legend.position = "top",
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plMutBenefit

# heatmaps

ch <- c(0.5, 1.5)
chm <- seq(1.5, 10, length.out = 150) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- seq(4, 12, length.out = 150)
dm <- 1

in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

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

plHeatMap <- ggplot(out_data,
       aes(x = chm, y = cm, fill = Outcome)) +
  geom_tile() + theme_classic() +
  scale_fill_manual(values = c("green", "orange", "red")) +
  labs(x = expression(atop("Added Host Colonization", "from the Mutalist" ~ (c[hm]))),
       y = expression("Mutualist Colonization" ~ (c[m])),
       fill = "") +
  facet_wrap(~ch, labeller = label_bquote(rows = c[h] == .(ch)), nrow = 2) +
  ggtitle("B") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.background = element_rect(fill = NA),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plHeatMap

jpeg("./figs/Fig3MutBenefit.jpeg",
     width = 3500, height = 2000, res = 300)
grid.arrange(plMutBenefit, plHeatMap,
            layout_matrix = matrix(c(1, 1, 1, 1, 2, 2, 2), nrow = 1, ncol = 7))
dev.off()




