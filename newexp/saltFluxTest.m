function saltFluxTest()
  bath = -410
  experiment_parameters = struct;
  experiment_parameters.tcline_atshelf_depth = bath-200;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.shelf_depth = 550;
  experiment_parameters.cavity_depth = -300;
  experiment_parameters.cavity_width = 150;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  experiment_parameters.rng_seed = 10;
  experiment_parameters.random_amplitude = 300;
  experiment_parameters.saltflux = false;
  experiment_parameters.rbcs_temp = false;

  currentFolder = pwd;
  experiment_parameters.tcline_atshelf_depth = bath+150;
  newexp('experiments/saltflux-explore/','above',experiment_parameters);
  cd('../experiments/saltflux-explore/above');
  system('sh upload_to_cluster.sh')
  cd(currentFolder);

end
