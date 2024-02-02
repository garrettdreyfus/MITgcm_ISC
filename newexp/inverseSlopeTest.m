function inverseSlopeTest()
  %% Depth at this seed is 465
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.shelf_depth = 700;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 32;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 250;
  experiment_parameters.saltflux = false;
  experiment_parameters.rbcs_temp = false;
  experiment_parameters.cavity_depth = 0;
  experiment_parameters.cavity_width = 150;
  experiment_parameters.yicefront = 150;
  currentFolder = pwd;
  bath= -495;
  indexs=[1];
  tcline_heights = [ -300, 0, 125 ]
  for i = indexs
    for k = tcline_heights
        experiment_parameters.tcline_atshelf_depth = bath+k;
        path_part1 = convertStringsToChars(strcat("experiments/inverse-GLIB-explore-",int2str(rng_seed),"/"));
        path_part2 = convertStringsToChars(strcat("d",int2str(experiment_parameters.cavity_depth),"at",int2str(k)));
        full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
        newexp(path_part1,path_part2,experiment_parameters);
        cd(full_path);
        system('sh upload_to_cluster.sh')
        cd(currentFolder);
    end
  end
end
