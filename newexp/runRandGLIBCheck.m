function runRandGLIBCheck()
  %% Depth at this seed is 465
  bath = -585;
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.shelf_depth = 650;
  experiment_parameters.cavity_depth = -300;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 16;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 150;
  experiment_parameters.saltflux = false;
  experiment_parameters.rbcs_temp = false;

  currentFolder = pwd;
  depths = [ -200 -125 -50 0 50 125 200] ;
  for k = depths
    experiment_parameters.tcline_atshelf_depth = bath+k;
    path_part1 = convertStringsToChars(strcat("experiments/GLIB-explore-",int2str(rng_seed),"/"));
    path_part2 = convertStringsToChars(strcat("at",int2str(k)));
    full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
    newexp(path_part1,path_part2,experiment_parameters);
    cd(full_path);
    system('sh upload_to_cluster.sh')
    cd(currentFolder);
  end
end