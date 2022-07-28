function runSlopeCheck()
  %% Depth at this seed is 465
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.shelf_depth = 650;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 18;
  experiment_parameters.cavity_width = 150;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 250;
  experiment_parameters.saltflux = false;
  experiment_parameters.rbcs_temp = false;
  depths = [-125 0 125] ;
  experiment_parameters.cavity_depth = -450;

  bath = -480
  experimen_parameters.bath = bath;
  currentFolder = pwd;
  for k = depths
    experiment_parameters.tcline_atshelf_depth = bath+k;
    path_part1 = convertStringsToChars(strcat("experiments/biggerslope-GLIB-explore-",int2str(rng_seed),"/"));
    path_part2 = convertStringsToChars(strcat("at",int2str(k)));
    full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    newexp(path_part1,path_part2,experiment_parameters);
    cd(full_path);
    system('sh upload_to_cluster.sh')
    cd(currentFolder);
  end
  currentFolder = pwd;
  depths = [-125 0 125] ;

  experiment_parameters.cavity_depth = -100;

  bath = -440
  experimen_parameters.bath = bath;
  for k = depths
    experiment_parameters.tcline_atshelf_depth = bath+k;
    path_part1 = convertStringsToChars(strcat("experiments/smallerslope-GLIB-explore-",int2str(rng_seed),"/"));
    path_part2 = convertStringsToChars(strcat("at",int2str(k)));
    full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    newexp(path_part1,path_part2,experiment_parameters);
    cd(full_path);
    system('sh upload_to_cluster.sh')
    cd(currentFolder);
  end
end
