def do_snapshot_espresso_system(espresso_system, filename):
      """
      Uses espresso visualizator for creating a snapshot of the current state of the espresso_system
      """ 
      from espressomd import visualization
      

      visualizer = visualization.openGLLive(
            espresso_system, bond_type_radius=[0.3], particle_coloring='type', draw_axis=False, background_color=[1, 1, 1],
      particle_type_colors=[[1.02,0.51,0], # Brown
                        [1,1,1],  # Grey
                        [2.55,0,0], # Red
                        [0,0,2.05],  # Blue
                        [0,0,2.05],  # Blue
                        [2.55,0,0], # Red
                        [2.05,1.02,0]]) # Orange
      visualizer.screenshot(filename)

      return
