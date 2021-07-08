import argparse
import ast
import configparser

from dms.tile import Tile

def bounded_number(converter, low=None, high=None):
    def f(s, converter=converter):
        try:
            x = converter(s)
        except:
            raise argparse.ArgumentTypeError(f"Failed to convert '{s}'")
        if (low is not None and x < low) or (high is not None and x > high):
            raise argparse.ArgumentTypeError(f'Invalid value: {x}')
        return x
    return f

non_negative_int = bounded_number(int, low=0)

def yes_or_no(s):
    s = s.lower()
    if s in ['true', 't', 'yes', 'y']:
        return True
    elif s in ['false', 'f', 'no', 'n']:
        return False
    else:
        raise argparse.ArgumentTypeError(
            f'Invalid value: {s}. Must be one of: True, T, Yes, Y,'
            ' False, F, No, N. (Not case sensitive.)')

def maybe_quoted_string(s):
    if len(s) > 0 and s[0] in ["'", '"']:
        return ast.literal_eval(s)
    else:
        return s

CLI_ONLY_ARGUMENTS = {
    'config' : {
        'help' : 'Configuration file to use.',
        'type' : str,
        'required' : True
    },
}

ARGUMENTS = {
    'use_multiprocessing' : {
        'help' : ('Should the multiprocessing module be used to process FASTQ'
                  ' files using multiple cores?'),
        'type' : yes_or_no,
        'default' : True
    },
    'fastq_file_dir' : {
        'help' : 'Directory where FASTQ files are located.',
        'type' : maybe_quoted_string,
        'default' : ''
    },
    'output_dir' : {
        'help' : 'Directory to place output files.',
        'type' : maybe_quoted_string,
        'default' : ''
    },
    'max_mismatches' : {
        'help' : ('Maximum number of mismatches a read can have'
                  ' without being discarded.'),
        'type' : non_negative_int,
        'default' : None
    },
    'min_quality' : {
        'help' : 'Reads with any quality score lower than this are discarded.',
        'type' : non_negative_int,
        'default' : None
    },
    'min_ref_counts' : {
        'help' : ('Variants with fewer than this number of reference counts'
                  ' will be discarded.'),
        'type' : non_negative_int,
        'default' : 5
    },
    'pseudocount' : {
        'help' : ('Pseudocount to use for mutations that were observed in the'
                  ' reference but not the selected population.'),
        'type' : non_negative_int,
        'default' : 1
    }
}

PARAMS_NAME = 'Parameters'
TILE_PREFIX = 'Tile:'
SAMPLES_NAME = 'Samples'
EXPERIMENTS_NAME = 'Experiments'
PROTEINS_NAME = 'Proteins'

def parse_params(arguments, config):
    params = argparse.Namespace()
    if not config.has_section(PARAMS_NAME):
        return params
    for param in config.options(PARAMS_NAME):
        if param not in arguments:
            raise ValueError(f'unknown option: {param}')
        try:
            raw = config.get(PARAMS_NAME, param)
            value = arguments[param]['type'](raw)
        except argparse.ArgumentTypeError as e:
            print(f'Invalid value for option {param}: {raw}')
            raise e
        setattr(params, param, value)
    return params

def parse_tile(config, section):
    kwargs = dict(wt_seq=ast.literal_eval(config.get(section, 'wt_seq')),
                  cds_start=config.getint(section, 'cds_start'),
                  cds_end=config.getint(section, 'cds_end'),
                  first_aa=config.getint(section, 'first_aa'))
    if config.has_option(section, 'positions'):
        kwargs['positions'] = ast.literal_eval(config.get(section, 'positions'))
    return Tile(**kwargs)

def parse_tiles(config):
    tiles = {}
    for section in config.sections():
        if section.startswith(TILE_PREFIX):
            tiles[section[len(TILE_PREFIX):]] = parse_tile(config, section)
    if len(tiles) == 0:
        raise argparse.ArgumentTypeError('config does not define any tiles.')
    return tiles

def parse_samples(config, tiles):
    if not config.has_section(SAMPLES_NAME):
        raise argparse.ArgumentTypeError('config does not have a [Samples]'
                                         ' section.')
    samples = {}
    for name, value in config.items(SAMPLES_NAME):
        elements = ast.literal_eval(value)
        if len(elements) != 3:
            raise ValueError(f'sample {name} does not specify a tile and two'
                             f' files.')
        tile, file1, file2 = elements
        if tile not in tiles:
            raise ValueError(f'sample {name} specifies undefined tile {tile}.')
        samples[name] = (tile, (file1, file2))
    return samples

def parse_experiments(config, samples):
    if not config.has_section(EXPERIMENTS_NAME):
        raise argparse.ArgumentTypeError('config does not have an [Experiments]'
                                         ' section.')
    experiments = {}
    for name, value in config.items(EXPERIMENTS_NAME):
        pair = ast.literal_eval(value)
        if len(pair) != 2:
            raise ValueError(f'experiment {name} does not specify two samples.')
        for sample in pair:
            if sample not in samples:
                raise ValueError(f'sample {sample} is not defined.')
        if samples[pair[0]][0] != samples[pair[1]][0]:
            raise ValueError(f'experiment {name} specifies samples with'
                             ' different tiles.')
        experiments[name] = tuple(pair)
    return experiments

def parse_proteins(config, tiles, samples, experiments):
    if not config.has_section(PROTEINS_NAME):
        raise argparse.ArgumentTypeError('config does not have a [Proteins]'
                                         ' section.')
    proteins = {}
    for name, value in config.items(PROTEINS_NAME):
        exps = ast.literal_eval(value)
        if len(exps) == 0:
            raise ValueError(f'protein {name} does not specify any'
                             ' experiments.')
        positions = {}
        for exp in exps:
            if exp not in experiments:
                raise ValueError(f'experiment {exp} is not defined.')
            sample_name = experiments[exp][0]
            tile_name = samples[sample_name][0]
            for pos in tiles[tile_name].positions:
                if pos in positions:
                    raise ValueError(f'protein {name}: includes position {pos}'
                                     f' from experiment {positions[pos]} and'
                                     f' experiment {exp}')
                else:
                    positions[pos] = exp
        proteins[name] = tuple(exps)
    if len(proteins) == 0:
        raise ValueError('config does not specify any proteins.')
    return proteins

def parse_config(path, args):
    """Parse a configuration file.

    Returns (params, tiles, samples, experiments). params is a dict
    mapping param_name -> value.  tiles is a dict mapping tile_name ->
    Tile. samples is a dict mapping sample_name -> (tile_name, (path1,
    path2)), where path1 and path2 are the forward and reverse
    paired-end reads for a single sample. experiments is a dict
    mapping experiment_name -> (ref_sample_name, sel_sample_name).
    """
    config = configparser.ConfigParser(strict=True)
    config.optionxform = str # make the parser case-sensitive
    config.read(path)
    params = parse_params(args, config)
    tiles = parse_tiles(config)
    samples = parse_samples(config, tiles)
    experiments = parse_experiments(config, samples)
    proteins = parse_proteins(config, tiles, samples, experiments)
    return params, tiles, samples, experiments, proteins

def make_arg_parser(cli_only_args, args):
    parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_args.items():
        parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    for name, kwargs in args.items():
        parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    return parser

def parse_args_and_read_config(argv,
                               cli_only_arguments=CLI_ONLY_ARGUMENTS,
                               arguments=ARGUMENTS):
    # First, parse CLI-only arguments.
    cli_only_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_arguments.items():
        cli_only_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    namespace, remaining_argv = cli_only_parser.parse_known_args(argv)
    # Next, get the configuration file and parse that.
    params, tiles, samples, experiments, proteins = \
        parse_config(namespace.config, arguments)
    # Finally, read everything else from the commandline, overriding
    # anything also specified in the config file Parameters block.
    arg_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in arguments.items():
        arg_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    arg_parser.parse_args(remaining_argv, params)
    return params, tiles, samples, experiments, proteins
