# call functions file
source("0_Functions.R")

# heatmaps for timescales

# initialize parameter values. Choose range of values for additional host colonization rate due to mutualist settlement and mutualist colonization rate. Choose 2 values for intrinsic host colonization rate. 
ch <- 0.25
chm <- seq(1.75, 20, length.out = 200) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

cp <- 5
dp <- 1

cm <- c(11, 12, 13)
dm <- 1

tau <- seq(1, 15, length.out = 200)
dm <- tau

# Find all possible combinations of parameters
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

in_pars$cp <- in_pars$cp * in_pars$dm
in_pars$dp <- in_pars$dp * in_pars$dm
in_pars$cm <- in_pars$cm * in_pars$dm

iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

# initialize output dataframe
out_data <- data.frame()

out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind') %dofuture% {
    
    # intiialize current iteration's parameter values
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)
    
    # solve for equilibria of each population
    root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
    h_soln <- root_soln$root
    mp_soln <- GetPM(h_soln, cur_pars)
    p_soln <- mp_soln[1]
    m_soln <- mp_soln[2]
    
    # extract leading eigenvalue from the Jacobian matrix to determine stability
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    cur_eig <- GetEig(J)
    
    # Determine feasibility and stability of equilbria
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    cur_stable <- cur_eig < 0
    
    # Label the feasibility and stability accordingly
    cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                          ifelse(cur_stable, "\nFeasible\nand stable\n",
                                 "\nFeasible\nbut unstable\n"))
    
    
    cur_params$unscaled_cm <- round(cur_params$cm / cur_params$dm, 3)
    cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
    cur_dyn
  }

# Weird stuff I had to do to change the labelling order to descend
out_data <- mutate(out_data, ch_label = paste0("c[h] == ", ch))

# create heatmap to visualize output data. X-axis represents additional host colonization due to mutualist settlement while the y-axis represents mutualist (colonizer) colonization rate
# Each panel indicates a different intrinsic host colonization rate
# Blue indices feasible and stable, yellow indicates feasible and unstable, and red indicates infeasible. 
plTimescales <- ggplot(out_data,
                    aes(x = chm, y = dm, fill = Outcome)) +
  geom_tile() + theme_classic() +
  scale_fill_manual("Coexistence\nstatus:", values = c("#0072B2", "#F0E442", "#D55E00")) +
  labs(x = expression(atop("Host Colonization due to the Mutualist" ~ (c[hm]))),
       y = expression("Timescale of bacterial to host dynamics" ~ (tau)),
       fill = "") +
   facet_wrap(~unscaled_cm, labeller = label_bquote(rows = c[m] == .(unscaled_cm)), nrow = 1) +
  ggtitle("A") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text = element_text( size = 10 ),
        panel.background = element_rect(fill = NA),
        panel.spacing = unit(1, "lines"),
        strip.background = element_blank(),
        text = element_text(size=15),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plTimescales

# analyzing fast bacterial dynamics

# Choose parameter values. Choose a range of values for the added colonization rate of host due to mutualist settlement. Choose three values for colonizer (mutualist) colonization rate
ch <- 0.25
chm <- seq(0.25+1e-10, 25, length.out = 200) # needs to be bigger than or equal to ch
dh <- 1
dhp <- 0

tau <- 100

cp <- 5 * tau
dp <- 1 * tau

cm <- 10 * tau
dm <- 1 * tau

# get all possible combinations of parameter values
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

pars <- list(ch = ch, chm = chm, dh = dh, dhp = dhp,
             cp = cp, dp = dp,
             cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choose initial condition standard deviation
ini_cond_sds <- c(0.001, 0.5)

# indiciate how long to run simulation and time step size
end_time <- 400
time_step <- 0.5

# initialize output data frame
out_data <- data.frame()

out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind') %:%
  foreach(
    ini_cond_sd = ini_cond_sds,
    .combine = 'rbind',
    .inorder = FALSE) %dofuture% {
      cur_params <- iterated_params[i,]
      cur_pars <- as.list(cur_params)
      
      # solve for equilibrium for each population
      root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
      h_soln <- root_soln$root
      mp_soln <- GetPM(h_soln, cur_pars)
      p_soln <- mp_soln[1]
      m_soln <- mp_soln[2]
      
      # determine feasibility
      cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
      
      # choose positive initial conditions to run simulation
      ini_state <- c(h_soln, p_soln, m_soln) + rnorm(n = 3, mean = 0, sd = ini_cond_sd)
      ini_state[ini_state <= 0] <- 0.0001
      names(ini_state) <- c("Host", "Pathogen", "Mutualist")
      
      # simulate the dynamics
      out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                                   end_time, time_step,
                                   fn = CoralMutPathDynamics)
      
      # extract last 100 time steps of the output
      out_dyn <- out_dyn %>%
        filter(time > (end_time - 100)) %>%
        melt(id.vars = c("time"))
      
      # extract leading eigenvalue for the Jacobian matrix to determine stability
      J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
      cur_eig <- GetEig(J)
      
      cur_params$IniCondSd <- ini_cond_sd
      
      # label outcomes accordingly with respective population, feasibility, and stability.
      cur_dyn <- cbind(cur_params, out_dyn) %>%
        mutate(pred_soln = case_when(variable == "Host" ~ h_soln,
                                     variable == "Mutualist" ~ m_soln,
                                     variable == "Pathogen" ~ p_soln)) %>%
        mutate(Eigenvalue = cur_eig,
               Stable = ifelse(cur_feas, ifelse(cur_eig < 0,
                                                "Feasible and stable",
                                                "Feasible but unstable"),
                               "Not feasible"))
      
      cur_dyn
    }

proc_data <- out_data %>%
  mutate(pred_soln = ifelse(pred_soln > 0, pred_soln, NA))

# Weird stuff I had to do to change the labelling order to descend
proc_data <- mutate(proc_data, cm_label = paste0("c[m] == ", cm))
proc_data$cm_label <- factor(proc_data$cm_label,
                             levels = c("c[m] == 15", "c[m] == 11", "c[m] == 10", "c[m] == 5"))

#ggplot(out_data,
#       aes(x = chm, y = Eigenvalue, color = as.factor(cm))) +
#  geom_point() + theme_classic() +
#  geom_hline(yintercept = 0) +
#  labs(x = "Added Host Colonization from the Mutualist",
#       y = "Eigenvalue",
#       color = "Mutualist\nColonization")

# plot output dataframe. X-axis represents additional host colonization rate from mutualist settlement while y-axis represents frequency of the respective population
# Blue indicates feasible and stable, yellow indicates feasible and unstable, and red indicates infeasible equilbria
# Each panel is a different population: host, mutualist (colonizer), or pathogen (competitor)
plFastBacteria <- ggplot(proc_data,
                       aes(x = chm, y = value)) +
  
  facet_grid(
    cols = vars(variable),
    labeller = label_parsed
  ) +
  geom_point(size = 3, alpha = 0.01, aes(color = as.factor(IniCondSd))) +
  theme_classic() +
  geom_line(linewidth = 2, aes(x = chm, y = pred_soln, color = Stable)) +
  scale_color_manual(breaks = c("0.001", "0.5", "Feasible and stable", "Feasible but unstable", "Not feasible"),
                     values = c("#000000", "#CC79A7", "#0072B2", "#F0E442", "#D55E00")) +
  guides(color = guide_legend(override.aes = list(linetype = 1, alpha = 1))) +
  labs(x = expression("Host Colonization due to the Mutualist" ~ (c[hm])),
       y = "Frequency",
       color = "Noise added to the Initial Conditions") +
  ggtitle("B") +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        legend.position = "top",
        strip.background = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plFastBacteria



# writing out graph
jpeg("./figs/Fig4Timescales.jpeg",
     width = 3500, height = 2900, res = 300)
grid.arrange(plTimescales, plFastBacteria, nrow = 2)
dev.off()





