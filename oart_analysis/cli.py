"""Default to bundled plotting data; opt into a local workbook with --source-data."""
import argparse
from pathlib import Path
import importlib.metadata
from .io import read_source, write_json
from .analysis import analyse
from .statistics import precision_sample_size


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('all','analyse','plot','metrics','sample-size'))
    parser.add_argument('--source-data',type=Path,default=None,
                        help='Use this local Source_data.xlsx instead of bundled plotting data.')
    parser.add_argument('--output',type=Path,default=Path('outputs'))
    parser.add_argument('--figure',choices=('5','6','S4','S5','S6','S7'))
    parser.add_argument('--config',type=Path,default=Path('config/imaging.local.json'))
    args=parser.parse_args(argv)
    if args.command=='sample-size':
        print(precision_sample_size())
        return
    if args.command=='metrics':
        from .imaging import evaluate_config
        count=evaluate_config(args.config,args.output)
        print(f'Evaluated {count} ROI pairs. Results: {args.output}')
        return
    if args.command=='plot' and not args.figure: parser.error('plot requires --figure')
    if args.command=='analyse' and args.source_data is None:
        parser.error('Full patient-level analysis requires --source-data /path/to/Source_data.xlsx. Bundled CSVs support figure reproduction; use all or plot without --source-data.')
    from .plot_data import from_workbook, load_bundled, validate
    if args.source_data is not None:
        if args.source_data.resolve().is_relative_to(args.output.resolve()):
            parser.error('Keep the input workbook outside the generated-output directory')
        print(f'Data source: explicitly selected workbook ({args.source_data.name}).',flush=True)
        source=read_source(args.source_data)
        report=analyse(source,args.output/'analysis')
        print(f"Verified {report['table_cells_verified']} table entries and {report['composite_values_verified']} patient OAR composites.")
        data, provenance=from_workbook(source)
    else:
        print('Data source: bundled plotting CSVs. No workbook is required or automatically searched.',flush=True)
        data, provenance=load_bundled()
        print('Plotting subset verified. Full patient-level retrospective analysis was not performed.',flush=True)
    provenance.update(command=args.command,figure=args.figure,validation=validate(data))
    write_json(provenance,args.output/'data_source.json')
    if args.command in ('all','plot'):
        from .plotting import render,FIGURES
        for figure in FIGURES if args.command=='all' else (args.figure,):
            print(f'Rendering Figure {figure}...',flush=True)
            folder=render(data,figure,args.output/'figures')
            write_json(provenance,folder/'data_source.json')
    versions={}
    for name in ('numpy','scipy','pandas','openpyxl','matplotlib','Pillow'):
        versions[name]=importlib.metadata.version(name)
    write_json(versions,args.output/'runtime_versions.json')
    print(f'Completed. Local outputs: {args.output}')


if __name__=='__main__':
    main()
