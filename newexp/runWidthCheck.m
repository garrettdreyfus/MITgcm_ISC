function runWidthCheck()
  %% Depth at this seed is 465
  experiment_parameters = struct;
  experiment_parameters.tcline_deltaz = 100;
  experiment_parameters.shelf_depth = 650;
  experiment_parameters.trough_depth = 0;
  experiment_parameters.rand_topo = true;
  experiment_parameters.monitor_freq = 6;
  rng_seed = 32;
  experiment_parameters.rng_seed = rng_seed;
  experiment_parameters.random_amplitude = 250;
  experiment_parameters.yicefront = 150;
  experiment_parameters.saltflux = false;
  experiment_parameters.rbcs_temp = false;
  experiment_parameters.cavity_depth = -300;
  currentFolder = pwd;
  cavity_widths = [ 50,100,250];
  baths= [-465 -470 -480];
  indexs=[1 2 3];
  tcline_heights = [ 125, 0, -300 ]
  for i = indexs
    w = cavity_widths(i);
    bath = baths(i);
    for k = tcline_heights
        experiment_parameters.tcline_atshelf_depth = bath+k;
        experiment_parameters.cavity_width = w;
        path_part1 = convertStringsToChars(strcat("experiments/widthexprevertnokpp-GLIB-explore-",int2str(rng_seed),"/"));
        path_part2 = convertStringsToChars(strcat("at",int2str(k),"w",int2str(w)));
        full_path = convertStringsToChars(strcat("../",path_part1,path_part2));
        newexp(path_part1,path_part2,experiment_parameters);
        cd(full_path);
        system('sh upload_to_cluster.sh')
        cd(currentFolder);
    end
  end
end
