def parse_args():
    from argparse import ArgumentParser
    import os

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('prmtop', type=os.path.realpath,
            help='Amber prmtop file')
    parser.add_argument('-c', '--coordinate-file', type=os.path.realpath,
            help='Amber coordinate file (rst7 or crd)')
    parser.add_argument('-o', '--output', type=os.path.realpath,
            required=True,
            help='Output file (dms, mae)')
    return parser.parse_args()

def main():
    import msys

    args = parse_args()
    mol = msys.LoadPrmTop(args.prmtop)
    if args.coordinate_file:
        msys.ReadCrdCoordinates(mol, args.coordinate_file)
    msys.Save(mol, args.output)

