theme_plot <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    text       = element_text(size = 16),
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16),
    
    panel.grid.major = element_line(color = scales::alpha("black", 0.2)),
    panel.grid.minor = element_line(color = scales::alpha("black", 0.1)),
    
    legend.position = "none"
  )
