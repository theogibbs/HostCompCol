# Call functions file
source("0_Functions.R")

# Feasibility and stability of coexistence as pathogen virulence and mutualist colonization rates are varied ----

ch <- c(2, 3, 4)
chm <-  ch # needs to be bigger than or equal to ch
dh <- 1
dhp <- seq(0, 2, length.out = 200)

cp <- 5
dp <- 1

cm <- seq(1e-10, 16, length.out = 200)
dm <- 1

# Get all possible combinations of parameter values
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)
in_pars$chm <- in_pars$ch

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))

# initialize output dataframe
out_data <- data.frame()
plan(multisession)
out_data <- foreach(
  i = 1:nrow(iterated_params),
  .combine = 'rbind',
  .inorder = FALSE) %dofuture% {

    # initialize current iteration's parameter values
    cur_params <- iterated_params[i,]
    cur_pars <- as.list(cur_params)

    # solve for equilibrium values of each population
    root_soln <- uniroot(f = PredEq, interval = c(0, 1), cur_pars)
    h_soln <- root_soln$root
    mp_soln <- GetPM(h_soln, cur_pars)
    p_soln <- mp_soln[1]
    m_soln <- mp_soln[2]

    # extract leading eigenvalue from Jacobian matrix to determine stability
    J <- BuildJacobian(h_soln, p_soln, m_soln, cur_pars)
    cur_eig <- GetEig(J)

    # determine feasibility and stability of the equilibrium
    cur_feas <- (h_soln > 0) * (p_soln > 0) * (m_soln > 0)
    cur_stable <- cur_eig < 0

    # Label feasibility and stability accordingly
    cur_outcome <- ifelse(!cur_feas, "\nNot feasible\n",
                          ifelse(cur_stable, "\nFeasible\nand stable\n",
                                 "\nFeasible\nbut unstable\n"))
    
    cur_dyn <- cbind(cur_params, data.frame(Outcome = cur_outcome))
    cur_dyn
  }

# create heatmap from output dataframe. X-axis represents pathogenicity value while y-axis represents colonizer (mutualist) colonization rate. 
# Red represents infeasible equilibrium while blue represents a feasible and stable equilibrium.
plHeatMapDeath <- ggplot(out_data,
                         aes(x = dhp, y = cm, fill = Outcome)) +
  facet_wrap(~ch, labeller = label_bquote(c[h] == .(ch))) +
  geom_tile() + theme_classic() +
  labs(x = expression("Added Host Mortality from the Pathogen" ~ (d[hp])),
       y = expression("Mutualist Colonization"~(c[m])),
       fill = "") +
  scale_fill_manual("Coexistence\nstatus:", values = c("#0072B2", "#D55E00")) +
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
plHeatMapDeath # Figure 2 B in main text

# Example trajectories as pathogen virulence is varied ----
ch <- 2
chm <- ch # needs to be bigger than or equal to ch
dh <- 1
dhp <- c(0, 0.25, 0.5)

cp <- 5
dp <- 1

cm <- 6
dm <- 1

# create all possible combinations of parameter values
in_pars <- crossing(ch = ch, chm = chm, dh = dh, dhp = dhp,
                    cp = cp, dp = dp,
                    cm = cm, dm = dm)

num_repl <- 1
iterated_params <- bind_rows(replicate(num_repl, in_pars, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# Indicate length of simulation and time step size
end_time <- 50
time_step <- 0.01

# initialize output dataframe
dyn_data <- data.frame()

for(i in 1:nrow(iterated_params)) {

  # initialize current iteration's parameter values
  cur_params <- iterated_params[i,]
  cur_pars <- as.list(cur_params)

  # set initial conditions for the simulation
  ini_state <- c(0.2, 0.05, 0.05)
  names(ini_state) <- c("Host", "Pathogen", "Mutualist")

  # simulate the dynamics 
  out_dyn <- IntegrateDynamics(ini_state, cur_pars,
                               end_time, time_step,
                               fn = CoralMutPathDynamics)
  
  cur_dyn <- cbind(cur_params, out_dyn)
  
  dyn_data <- rbind(dyn_data, cur_dyn)
  
}

# Extract the pathogenicity rate and densities for each population from the output dataframe
melt_dyn <- dyn_data %>%
  select(c("dhp", "time", "Host", "Pathogen", "Mutualist")) %>%
  melt(id.vars = c("time", "dhp"))

# plot the dynamics. Time on the x-axis and proportion of space occupied on the y-axis
plDyn <- ggplot(melt_dyn, aes(x = time, y = value, color = variable)) +
  geom_line(linewidth = 1) + theme_classic() +
  facet_wrap(~dhp, scales = "free", labeller = label_bquote(cols = d[hp] == .(dhp))) +
  labs(x = "Time", y = "Frequency", color = "") +
  ggtitle("A") +
  scale_color_viridis_d() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.05)) +
  theme(text = element_text(size=15),
        legend.text=element_text(size = 15),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.title.position = "plot",
        plot.caption.position =  "plot")
plDyn # Figure 2 A in the main text

# Save figure ----
jpeg("./figs/Fig2PathDeath.jpeg",
     width = 3000, height = 2000, res = 300)
grid.arrange(plDyn, plHeatMapDeath, nrow = 2)
dev.off()
